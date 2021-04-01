#!/usr/bin/py
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
    shell: "cutadapt -m 15 -a {params.adapter1} -A {params.adapter2} -n 2 -o {output.trimmed_read1} -p {output.trimmed_read2} {input.all_read1} {input.all_read2} >{output.trimmed_stats}"

rule map:
    input:
        p1 = rules.trim.output.trimmed_read1,
        p2 = rules.trim.output.trimmed_read2
    output:
        mapped_bam_file = config["dir_names"]["mapped_dir"] + "/{sample_id}.bam",
    version: config["tool_version"]["bwa"]
    params:
        threads = config["params"]["bwa"]["threads"],
        map_all = config["params"]["bwa"]["all"],
        reference = config["params"]["bwa"]["bwa_reference"]
    shell:
        """
        bwa mem -t {params.threads} {params.reference} {input.p1} {input.p2} | samtools view -F 4 -Sb |  samtools sort -T {output.mapped_bam_file}.align -o {output.mapped_bam_file}
        """

rule ivar_filter:
    input:
        sorted_bam = rules.map.output.mapped_bam_file,
        primer_bed = config["params"]["ivar"]["primer_bed"]
    output:
        filtered_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.filtered.bam"
    shell:
        """
        ivar trim -i {input.sorted_bam} -b {input.primer_bed} -p {output.filtered_bam_file}
        
        ### IF NEB KIT UNCOMMENT BELOW and comment ABOVE LINE
        ### ivar trim -e -k -i {input.sorted_bam} -b {input.primer_bed} -p {output.filtered_bam_file}
        """

rule second_sort:
    input:
        filtered_bam = rules.ivar_filter.output.filtered_bam_file
    output:
        second_sorted_bam_file = config["dir_names"]["filtered_dir"] + "/{sample_id}.filtered.sorted.bam"
    shell:
        """
        samtools view -hb {input.filtered_bam} | samtools sort -T {input.filtered_bam}.tmp -o {output.second_sorted_bam_file}
        """

rule pileup:
    input:
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        reference = config["params"]["bwa"]["bwa_reference"]
    output:
        pileup = config["dir_names"]["mpileup_dir"] + "/{sample_id}.mpileup"
    params:
        depth = config["params"]["mpileup"]["depth"],
        min_base_qual = config["params"]["varscan"]["snp_qual_threshold"]
    shell:
        """
        bcftools mpileup -A \
            -Q {params.min_base_qual} \
            --max-depth 5000000 -L 5000000\
            -f {input.reference} \
            {input.second_sorted_bam} > {output.pileup} 
        """

rule call_snps:
    input:
        pileup = rules.pileup.output.pileup,
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        reference = config["params"]["bwa"]["bwa_reference"]
    output:
        vcf = config["dir_names"]["varscan_dir"] + "/{sample_id}.vcf.gz"
    params:
        depth = config["params"]["mpileup"]["depth"],
        min_base_qual = config["params"]["varscan"]["snp_qual_threshold"],
        min_base_cov = config["params"]["varscan"]["min_cov"]
    shell:
        """
        bcftools mpileup -A \
            -a "INFO/AD,INFO/ADF,INFO/ADR,FORMAT/ADF,FORMAT/ADR,FORMAT/SP"\
            -Q {params.min_base_qual} \
            -f {input.reference} -L 5000000 --max-depth 5000000 \
            -Ou \
            {input.second_sorted_bam} |\
         bcftools call -Ou -mv |\
         bcftools norm -f {input.reference} -Ou |
         bcftools filter --include '(TYPE="INDEL" && IMF >.3 && IDV > 30) || (TYPE="SNP" && DP > {params.min_base_cov})' -Oz -o {output.vcf}
        """

rule call_consensus_snps:
    input:
        vcf = rules.call_snps.output.vcf
    output:
        consensus_vcf = config["dir_names"]["varscan_dir"] + "/{sample_id}.consensus.vcf.gz"
    shell:
        """
        bcftools filter --exclude '(AD[0])/ (AD[0] + AD[1]) >= 0.5' {input.vcf} -Oz -o {output.consensus_vcf}
        """

rule index_bcf:
    input:
        vcf = rules.call_snps.output.vcf,
        consensus_vcf = rules.call_consensus_snps.output.consensus_vcf
    output:
        index = config["dir_names"]["varscan_dir"]+"/{sample_id}.vcf.gz.tbi",
        consensus_index = config["dir_names"]["varscan_dir"]+"/{sample_id}.consensus.vcf.gz.tbi"
    shell:
        """
        tabix -f {input.vcf} > {output.index};
        tabix -f {input.consensus_vcf} > {output.consensus_index};
        """

rule bedtools_mask:
    input:
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        vcf = rules.call_snps.output.vcf
    output:
        bed_file = config["dir_names"]["varscan_dir"]+"/{sample_id}.bed"
    shell:
        """
        bedtools genomecov -ibam {input.second_sorted_bam} -bga |\
        awk '{{if($4 < 10) print $_}}' |\
        bedtools intersect -v -a - -b {input.vcf} > {output.bed_file}
        """

rule mask_consensus:
    input:
        vcf = rules.call_consensus_snps.output.consensus_vcf,
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        index = rules.index_bcf.output.consensus_index,
        bed_file = rules.bedtools_mask.output.bed_file
    output:
        consensus_genome = config["dir_names"]["consensus_dir"]+"/{sample_id}.masked_consensus.fasta"
    params:
        reference = config["params"]["bwa"]["bwa_reference"]
    shell:
        """
        bcftools consensus -f {params.reference} -m {input.bed_file} {input.vcf} > {output.consensus_genome}
        """
