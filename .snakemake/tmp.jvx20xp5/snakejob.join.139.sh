#!/bin/sh
# properties = {"type": "single", "rule": "join", "local": false, "input": ["outputs/unzip/COIII_3B_MRB7260_plustet1_R1.fq", "outputs/unzip/COIII_3B_MRB7260_plustet1_R2.fq"], "output": ["outputs/join/COIII_3B_MRB7260_plustet1.assembled.fastq"], "wildcards": {"sample_id": "COIII_3B_MRB7260_plustet1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 139, "cluster": {"cluster": "faculty", "partition": "gbc", "qos": "gbc", "time": "02:00:00", "nodes": 1, "ntasks-per-node": 1}}
cd /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv && \
/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/bin/python3 \
-m snakemake outputs/join/COIII_3B_MRB7260_plustet1.assembled.fastq --snakefile /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.jvx20xp5 outputs/unzip/COIII_3B_MRB7260_plustet1_R1.fq outputs/unzip/COIII_3B_MRB7260_plustet1_R2.fq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules join --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.jvx20xp5/139.jobfinished" || (touch "/projects/academic/lread/core-sequencing-runs/treat-scripting/development/treat-virtualenv/.snakemake/tmp.jvx20xp5/139.jobfailed"; exit 1)

