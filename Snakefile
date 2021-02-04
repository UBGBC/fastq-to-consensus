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
            -Q {params.min_base_qual} \
            -f {input.reference} -L 5000000 --max-depth 5000000 \
            -Ou \
            {input.second_sorted_bam} |\
         bcftools call -Ou -mv |\
         bcftools norm -f {input.reference} -Ou |
         bcftools filter --include '(TYPE="INDEL" && IMF >.3 && IDV > 30) || (TYPE="SNP" && DP > {params.min_base_cov})' -Oz -o {output.vcf}
        """

#rule call_snps_varscan:
#    input:
#        pileup = rules.pileup.output.pileup
#    output:
#        vcf = config["dir_names"]["varscan_dir"] +"/{sample_id}.varscan.vcf"
#    params:
#        min_cov = config["params"]["varscan"]["min_cov"],
#        snp_qual_threshold = config["params"]["varscan"]["snp_qual_threshold"],
#        snp_frequency = config["params"]["varscan"]["snp_frequency"]
#    shell:
#        """
#        java -jar tools/varscan/VarScan.v2.3.9.jar mpileup2cns \
#            {input.pileup} \
#            --min-coverage {params.min_cov} \
#            --min-avg-qual {params.snp_qual_threshold} \
#            --min-var-freq {params.snp_frequency} \
#            --strand-filter 0 \
#            --variants \
#            --output-vcf 1 > {output.vcf}
#        """

rule zip_vcf:
    input:
        vcf = rules.call_snps.output.vcf
    output:
        bcf = config["dir_names"]["varscan_dir"]+"/{sample_id}.vcf.gz"
    shell:
        """
        bgzip -c {input.vcf} > {output.bcf}
        """

rule index_bcf:
    input:
        vcf = rules.call_snps.output.vcf
    output:
        index = config["dir_names"]["varscan_dir"]+"/{sample_id}.vcf.gz.tbi"
    shell:
        """
        tabix {input.vcf}
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
        vcf = rules.call_snps.output.vcf,
        index = rules.index_bcf.output.index,
        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file,
        bed_file = rules.bedtools_mask.output.bed_file
    output:
        consensus_genome = config["dir_names"]["consensus_dir"]+"/{sample_id}.masked_consensus.fasta"
    params:
        reference = config["params"]["bwa"]["bwa_reference"]
    shell:
        """
        bcftools consensus -f {params.reference} -m {input.bed_file} {input.vcf} > {output.consensus_genome}
        """

#rule vcf_to_consensus:
#    input:
#        bcf = rules.zip_vcf.output.bcf,
#        index = rules.index_bcf.output.index,
#        ref = config["params"]["bwa"]["bwa_reference"]
#    output:
#        consensus_genome = config["dir_names"]["consensus_dir"]+"/{sample_id}.consensus.fasta"
#    shell:
#        """
#        cat {input.ref} | \
#            bcftools consensus {input.bcf} > \
#            {output.consensus_genome}
##        """

#rule create_bed_file:
#    input:
#        pileup = rules.pileup_to_consensus.output.pileup,
#        second_sorted_bam = rules.second_sort.output.second_sorted_bam_file
#    output:
#        bed_file = config["dir_names"]["varscan_dir"]+"/{sample_id}.bed"
#    params:
#        min_cov = config["params"]["varscan"]["min_cov"],
#        min_freq = config["params"]["varscan"]["snp_frequency"],
#        ref = config["params"]["bwa"]["bwa_reference"]
#    shell:
#        """
#        #samtools mpileup -a -d 5000000 -f {params.ref} {input.second_sorted_bam} > {input.pileup}.samtools;
#        python tools/seattleflu-scripts/create_bed_file_for_masking.py \
#            --pileup {input.pileup} \
#            --min-cov {params.min_cov} \
#            --min-freq {params.min_freq} \
#            --bed-file {output.bed_file}
#        """

#rule mask_consensus:
#    input:
#        consensus_genome = rules.pileup_to_consensus.output.consensus_genome,
#        low_coverage = rules.create_bed_file.output.bed_file
#    output:
#        masked_consensus = config["dir_names"]["consensus_dir"]+"/{sample_id}.masked_consensus.fasta"
#    shell:
#        """
#        bedtools maskfasta \
#            -fi {input.consensus_genome} \
#            -bed {input.low_coverage} \
#            -fo {output.masked_consensus}
#        """
