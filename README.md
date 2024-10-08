This pipeline is designed to process COVID-19 sequencing data generated at the UB Genomics and Bioinformatics Core Facility (http://ubnextgencore.buffalo.edu). We are operating on illumina MiSeq and NovaSeq short-read sequencing data. Modules loaded are custom to University at Buffalo Center for Computational Research (CCR) though are all open source and freely available. 

authors: Jonathan Bard (jbard at buffalo dot edu) , Brandon Marzullo , Dr. Natalie Lamb

<h2>Implementation</h2>
Snakemake drives this pipeline and submits individual processing jobs to our SLURM analysis cluster. Adapters are removed using Cutadapt, prior to mapping to the COVID-19 reference genome using bwa mem. Ivar then filteres our alignment files to mask ARCTIC primer locations, to avoid calling variants originating from primer oligonucleotides. Variant calling is handled bcftools with parameters -Q {quality filter} , -A (orphan reads), -L {max.depth for indel calling} and --max-depth {max.depth}, to ensure full use of our data. In addition, INDELs require the vcf flag IMF > .3 and  IDV > {min_base_cov} , while SNPs require DP > {min_base_cov}. These ensure we detect snp and INDELs with at least a specified minimum read depth. Bcftools commands: mpileup -> call -> norm -> filter. Lastly prior to consensus calling, a per-base-pair depth calculation using bedtools genomecov produces a low-depth masking file which is provided to bcftools consensus along with the reference, variant calls to ensure accurate consensus genome calls.

<h1> # POST-EASYBUILD CCR Directions </h1>

`salloc --qos general-compute  --nodes 1 --cpus-per-task 12 --mem 128G`

`srun --pty /bin/bash --login`

`module --initial_load restore`

`module use /projects/academic/gbcstaff/utils/modules/`

`module load gcc/11.2.0 openmpi/4.1.1 snakemake/7.18.2 gbc-anaconda/3.6 samtools/1.16.1 bowtie2/2.4.4 gbc-cutadapt/3.5 htslib/1.14 bcftools/1.14 bedtools/2.30.0 bwa/0.7.17`

`export LD_LIBRARY_PATH=/projects/academic/gbcstaff/utils/ivar-2024/htslib/lib/:/cvmfs/soft.ccr.buffalo.edu/versions/2022.05/easybuild/software/Core/gcccore/11.2.0/lib64/:/cvmfs/soft.ccr.buffalo.edu/versions/2022.05/easybuild/software/Core/gcccore/11.2.0/lib64/:/cvmfs/soft.ccr.buffalo.edu/versions/2023.01/compat/lib64/libc.so.6`

`export PATH=$PATH:/projects/academic/gbcstaff/utils/ivar-2024/bin`

`git clone https://github.com/UBGBC/fastq-to-consensus`

`source bin/activate`

`snakemake --latency-wait 120 -p -j 100 --profile slurm`


<h1># fastq-to-consensus step-by-step</h1>
Handles the preprocessing from illumina fastq files throught to consensus genome fasta as compared to a reference genome

<h3>Currently Loaded Modules:</h3>

  `module load gcc/11.2.0 gbc-ivar/1.2.4 gbc-samtools/1.10 gbc-bowtie2/2.4.1 gbc-cutadapt/3.5 htslib/1.2.1 gbc-bcftools/1.9 gbc-bedtools/2.29.1 bwa/0.7.17`
  
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

    `snakemake --latency-wait 120 -p -j 100 --profile slurm`  

7. Pipeline should result in a consensus.fasta file per sample



### NO LONGER NEEDED ### 
Consensus file fasta headers need to be renamed  --

`for i in $(ls *masked_consensus.fasta | cut -f1 -d"."); do echo "sed -i 's/>MN908947\.3/>hCov-19\/USA\/NY-$i\/2021/' $i.masked_consensus.fasta ";done > cmds.txt`


## Potential errors:
This is caused by duplicate samples in metadata sheet.
`TypeError: join() argument must be str or bytes, not 'Series'
Wildcards:`
