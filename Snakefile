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
    input:expand("{dir}/{sample_id}.sorted.bam", dir=config["dir_names"]["sorted_dir"],sample_id=sample_ids)
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
    version: config["tool_version"]["bowtie2"]
    params:
        threads = config["params"]["bowtie2"]["threads"],
        map_all = config["params"]["bowtie2"]["all"],
        reference = config["params"]["bowtie2"]["bowtie2_reference"]
    shell:
        """
        bowtie2 \
            -x {params.reference} \
            -1 {input.p1} \
            -2 {input.p2} \
            -P {params.threads} \
            {params.map_all} \
            --local 2> {output.bt2_log} | \
                samtools view -bSF4 - > {output.mapped_bam_file}
        """

rule sort:
    input:
        mapped_bam = rules.map.output.mapped_bam_file
    output:
        sorted_bam_file = config["dir_names"]["sorted_dir"] + "/{sample_id}.sorted.bam"
    shell:
        """
        samtools view \
            -bS {input.mapped_bam} | \
            samtools sort | \
            samtools view -h > {output.sorted_bam_file}
        """


