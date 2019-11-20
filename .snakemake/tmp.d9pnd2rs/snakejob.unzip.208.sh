#!/bin/sh
# properties = {"type": "single", "rule": "unzip", "local": false, "input": ["../C9K8C/fastq/CR3_REH1_plustet2_S32_L001_R1_001.fastq.gz", "../C9K8C/fastq/CR3_REH1_plustet2_S32_L001_R2_001.fastq.gz", "mkdir.done"], "output": ["outputs/unzip/CR3_REH1_plustet2_R1.fq", "outputs/unzip/CR3_REH1_plustet2_R2.fq"], "wildcards": {"sample_id": "CR3_REH1_plustet2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 208, "cluster": {"cluster": "faculty", "partition": "gbc", "qos": "gbc", "time": "00:30:00", "nodes": 1, "ntasks-per-node": 1}}
cd /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv && \
/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/bin/python3 \
-m snakemake outputs/unzip/CR3_REH1_plustet2_R1.fq --snakefile /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.d9pnd2rs ../C9K8C/fastq/CR3_REH1_plustet2_S32_L001_R1_001.fastq.gz ../C9K8C/fastq/CR3_REH1_plustet2_S32_L001_R2_001.fastq.gz mkdir.done --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules unzip --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.d9pnd2rs/208.jobfinished" || (touch "/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.d9pnd2rs/208.jobfailed"; exit 1)

