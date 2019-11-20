#!/bin/sh
# properties = {"type": "single", "rule": "unzip", "local": false, "input": ["../C9K8C/fastq/CR4_3B_MRB7260_mintet1_S9_L001_R1_001.fastq.gz", "../C9K8C/fastq/CR4_3B_MRB7260_mintet1_S9_L001_R2_001.fastq.gz", "mkdir.done"], "output": ["outputs/unzip/CR4_3B_MRB7260_mintet1_R1.fq", "outputs/unzip/CR4_3B_MRB7260_mintet1_R2.fq"], "wildcards": {"sample_id": "CR4_3B_MRB7260_mintet1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 177, "cluster": {"cluster": "faculty", "partition": "gbc", "qos": "gbc", "time": "00:30:00", "nodes": 1, "ntasks-per-node": 1}}
cd /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv && \
/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/bin/python3 \
-m snakemake outputs/unzip/CR4_3B_MRB7260_mintet1_R1.fq --snakefile /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.kw3n28f2 ../C9K8C/fastq/CR4_3B_MRB7260_mintet1_S9_L001_R1_001.fastq.gz ../C9K8C/fastq/CR4_3B_MRB7260_mintet1_S9_L001_R2_001.fastq.gz mkdir.done --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules unzip --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.kw3n28f2/177.jobfinished" || (touch "/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.kw3n28f2/177.jobfailed"; exit 1)

