<h1># fastq-to-consensus</h1>
Handles the preprocessing from illumina fastq files throught to consensus genome fasta as compared to a reference genome

Currently Loaded Modules:
  1) gbc-samtools/1.7   2) gbc-bowtie2/2.4.1

1. Navigate to the new flowcell data output.

2. git clone this repository 

    `git clone https://github.com/UBGBC/fastq-to-consensus`

3. Activate the python anaconda environment (testing on CCR 11-21-19)

    `source fastq-to-consensus/bin/activate` 

4. Edit the config.json file and cluster.json files


5. Ensure meta-data table contains all of the necessairy fields

** NOTE EXACT HEADERS HAVE TO BE ENFORCED or key errors will be thrown during processing**


6. Launch jobs

  The pipeline will utilize CCR resource to parallel execution.
  OTU table and statisics about merge rate, filter rate, hit rate wiil be placed under _table_

### The use of --latency-wait allows for SLURM to catch up writing the files and posting the file handles so Snakemake can see them.

    `snakemake --latency-wait 120 -p -j 100 --cluster-config cluster.json --cluster "sbatch --partition gbc --cluster faculty --qos gbc --account gbcstaff"`
   

7. Pipeline should result in a consensus.fasta file per sample
