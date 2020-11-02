This pipeline is a heavily modified aggregation of two github repositories:

https://github.com/seattleflu/
https://github.com/jebard/fastq-to-treat



<h1># fastq-to-consensus</h1>
Handles the preprocessing from illumina fastq files throught to consensus genome fasta as compared to a reference genome

<h3>Currently Loaded Modules:</h3>

  1) gcc/10.2.0   2) gbc-ivar/1.2.4   3) gbc-samtools/1.10   4) gbc-bowtie2/2.4.1   5) gbc-cutadapt/1.16   6) htslib/1.2.1   7) gbc-bcftools/1.9   8) gbc-bedtools/2.29.1   9) bwa/0.7.17
  
  module load gbc-ivar gbc-samtools gbc-bowtie2 gbc-cutadapt htslib gbc-bcftools gbc-bedtools bwa
  
  
<h3> Step-by-step of install and analysis </h3>
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
