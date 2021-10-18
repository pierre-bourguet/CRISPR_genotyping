# CRISPR_genotyping
A tool to genotype indels generated by multiplexed CRISPR-Cas9 using amplicon-sequencing. The script scans for indels up to 20bp in size and reports those found in at least 10% of all sequenced reads (frequency > 10%). Substitutions are ignored.

## Usage
sbatch main.sh -q fastq_directory -f fasta_file.fa -n study_name  
### Options
-q directory to your raw sequencing data in fq.gz format. If your file name has a lot of useless stuff in it, it will be kept in this pipeline. For clarity, use something like sampleID.fq.gz  
-f fasta file that contains your target genes (relative or absolute path).  
-n study name, a prefix that will be added to files and folders.  
The command should be run from the directory where the script files are located. The main program is to be run on the cluster.  
Maximum tested indel size and threshold for reported indels can be customized easily by changing argument values at the beginning of *main.sh*  
The interactive version (slower but not dependent on the slurm queue) is meant to be run in interactive mode (not recommended). This code can be useful if you want to run this locally though.  

## Output files
A directory with your study name is created in the current directory and contains all output files:  
* **total_frequency_of_indels.tsv** : for each sample and gene, adds up the frequencies of all indels that passed the frequency threshold (indel found in > 10% of the reads by default). Table is sorted from sample with highest sum of all frequencies to lowest. This table is used to select samples with the lowest frequency of WT alleles, but does not take into account whether indels are fixed (single indel with high frequency, typically close to 50% or 100%) or if they are a mixture of many different low-frequency mutations. I would recommend using this table to select lines rather than the file below.
* **frequency_of_top_coverage_indels.tsv** : for each sample and gene, selects the indel with the highest **coverage** (we would rather like highest frequency here, but in our first dataset we often found that highest frequency mutations are low coverage mutations which are likely noise) but reports its **frequency**. In case of indels with equal coverage, selects the one with highest frequency. Table is sorted from sample with highest sum of all indel frequencies to lowest (sum of the row values).
* **idxstats.tsv** Aggregated output from samtools idxstats. Reports the gene size in 1st column followed by number of reads at each gene for each sample. Used for quality control, coverage should be relatively homogeneous.
* **bam**: folder with the aligned reads (.bam) and index files (.bai). Needed to vizualize reads in IGV or another genome browser
* **mutations**: folder with all indels above the frequency threshold (10% by default) for each gene. Columns are: genename, position, reference sequence, coverage, indel, frequency.
* bwa_logs, indexes, mpileup folders : intermediate folders you can ignore

## Analysis
Workflow is as follow: use .tsv tables to find samples with the most indels. Check the most interesting samples in IGV by dragging bam files. The reason to check is that the script only tests indels up to 20bp, so substitutions and larger indels are missed (not so rare). IGV also helps to vizualize alleles with complex mutations.  
Small tip: for indels reported in the mutations folder, the reference sequence is copied from your input fasta file. You can quickly check if indels occur where they should by using lowercase for your CRISPR target sequence in the fasta file, while the rest of the sequence is in uppercase. If indels are detected in this target sequence, the reference base in "mutations" files will be lowercase.
