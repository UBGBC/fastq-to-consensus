import pandas as pd
import os

configfile: "config.json"
localrules: all, mkdir

df = pd.read_csv(config["meta_file"], sep='\t', header=0, index_col=0)
sample_ids = list(df.index)
df.index = sample_ids

def get_pair_gz(sample_id):
    dir = config["raw_fastq_gz_dir"]
    return tuple(os.path.join(dir, df.loc[str(sample_id), x]) for x in ('ForwardFastqGZ', 'ReverseFastqGZ'))

def get_forward_primer(sample_id):
    return df.loc[sample_id]["Adapter_1"]

def get_reverse_primer(sample_id):
    return df.loc[sample_id]["Adapter_2"]

rule all:
    input:expand("{dir}/{sample_id}.masked_consensus.fasta", dir=config["dir_names"]["consensus_dir"],sample_id=sample_ids)
    run:
        for sample in sample_ids:
            print("Wrapping up pipeline")
                   
rule mkdir:
    output: touch(config["file_names"]["mkdir_done"])
    params: dirs = list(config["dir_names"].values())
    shell: "mkdir -p {params.dirs}"

rule trim:
    input: 
        rules.mkdir.output,
        all_read1 = lambda wildcards: get_pair_gz(wildcards.sample_id)[0],
        all_read2 = lambda wildcards: get_pair_gz(wildcards.sample_id)[1]
    output: 
        trimmed_read1 = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R1.fastq.gz",
        trimmed_read2 = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.R2.fastq.gz",
        trimmed_stats = config["dir_names"]["trimmed_dir"] + "/{sample_id}.trimmed.stats"

    version: config["tool_version"]["cutadapt"]
    params:
        adapter1=lambda wildcards: get_forward_primer(wildcards.sample_id),
        adapter2=lambda wildcards: get_reverse_primer(wildcards.sample_id)
    shell: "cutadapt -m 15 -a {params.adapter1} -A {params.adapter2} -n 2 -o {output.trimmed_read1} -p {output.trimmed_read2} {input.all_read1} {input.all_read2} 2>{output.trimmed_stats}"

rule map:
    input:
        p1 = rules.trim.output.trimmed_read1,
        p2 = rules.trim.output.trimmed_read2
    output:
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam",
        bt2_log = config["dir_names"]["mapped_dir"] + "/{sample_id}.log"
    version: config["tool_version"]["bwa"]
    params:
        threads = config["params"]["bwa"]["threads"],
        map_all = config["params"]["bwa"]["all"],
        reference = config["params"]["bwa"]["bwa_reference"]
    shell:
        """
        bwa mem -t {params.threads} {params.reference} {input.p1} {input.p2} | samtools sort | samtools view -F 4 -o {output.mapped_bam_file}
        """

rule sort:
    input:
        mapped_bam = rules.map.output.mapped_bam_file
    output:
        sorted_bam_file = config["dir_names"]["sorted_dir"] + "/{sample_id}.sorted.bam"
    shell:
        """
        samtools view \
            -hS {input.mapped_bam} | \
            samtools sort | \
            samtools view -hb > {output.sorted_bam_file}
        """

rule ivar_filter:
    input:
        sorted_bam = rules.sort.output.sorted_bam_file,
        primer_bed = config["params"]["ivar"]["primer_bed"]
    output:
        filtered_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.filtered.bam"
    shell:
        """
        ivar trim -i {input.sorted_bam} -b {input.primer_bed} -p {output.filtered_bam_file}
        """

rule second_sort:
    input:
        filtered_bam = rules.ivar_filter.output.filtered_bam_file
    output:
        second_sorted_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.filtered.sorted.bam"
    shell:
        """
        samtools view \
            -hS {input.filtered_bam} | \
            samtools sort | \
            samtools view -hb > {output.second_sorted_bam_file}
        """

rule pileup:
    input:
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        reference = config["params"]["bowtie2"]["bowtie2_reference"]+".fasta"
    output:
        pileup = config["dir_names"]["mpileup_dir"] + "/{sample_id}.mpileup"
    params:
        depth = config["params"]["mpileup"]["depth"],
        min_base_qual = config["params"]["varscan"]["snp_qual_threshold"]
    shell:
        """
        samtools mpileup -a -A \
            -Q {params.min_base_qual} \
            -d {params.depth} \
            {input.second_sorted_bam} > {output.pileup} \
            -f {input.reference}
        """

rule call_snps:
    input:
        pileup = rules.pileup.output.pileup
    output:
        vcf = config["dir_names"]["varscan_dir"] +"/{sample_id}.vcf"
    params:
        min_cov = config["params"]["varscan"]["min_cov"],
        snp_qual_threshold = config["params"]["varscan"]["snp_qual_threshold"],
        snp_frequency = config["params"]["varscan"]["snp_frequency"]
    shell:
        """
        java -jar tools/varscan/VarScan.v2.3.9.jar mpileup2snp \
            {input.pileup} \
            --min-coverage {params.min_cov} \
            --min-avg-qual {params.snp_qual_threshold} \
            --min-var-freq {params.snp_frequency} \
            --strand-filter 1 \
            --output-vcf 1 > {output.vcf}
        """

rule zip_vcf:
    input:
        vcf = rules.call_snps.output.vcf
    output:
        bcf = config["dir_names"]["varscan_dir"]+"/{sample_id}.vcf.gz"
    shell:
        """
        bgzip {input.vcf}
        """

rule index_bcf:
    input:
        bcf = rules.zip_vcf.output.bcf
    output:
        index = config["dir_names"]["varscan_dir"]+"/{sample_id}.vcf.gz.csi"
    shell:
        """
        bcftools index {input}
        """

rule vcf_to_consensus:
    input:
        bcf = rules.zip_vcf.output.bcf,
        index = rules.index_bcf.output.index,
        ref = config["params"]["bowtie2"]["bowtie2_reference"]+".fasta"
    output:
        consensus_genome = config["dir_names"]["consensus_dir"]+"/{sample_id}.consensus.fasta"
    shell:
        """
        cat {input.ref} | \
            bcftools consensus {input.bcf} > \
            {output.consensus_genome}
        """

rule create_bed_file:
    input:
        pileup = rules.pileup.output.pileup
    output:
        bed_file = config["dir_names"]["varscan_dir"]+"/{sample_id}.bed"
    params:
        min_cov = config["params"]["varscan"]["min_cov"],
        min_freq = config["params"]["varscan"]["snp_frequency"]
    shell:
        """
        python tools/seattleflu-scripts/create_bed_file_for_masking.py \
            --pileup {input.pileup} \
            --min-cov {params.min_cov} \
            --min-freq {params.min_freq} \
            --bed-file {output.bed_file}
        """

rule mask_consensus:
    input:
        consensus_genome = rules.vcf_to_consensus.output.consensus_genome,
        low_coverage = rules.create_bed_file.output.bed_file
    output:
        masked_consensus = config["dir_names"]["consensus_dir"]+"/{sample_id}.masked_consensus.fasta"
    shell:
        """
        bedtools maskfasta \
            -fi {input.consensus_genome} \
            -bed {input.low_coverage} \
            -fo {output.masked_consensus}
        """
