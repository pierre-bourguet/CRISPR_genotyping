#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-00:05:00
#SBATCH --qos=short
#SBATCH --array=1-131
#SBATCH --output=logs/%A_%a_%x.out
#SBATCH --error=logs/%A_%a_%x.err
# above is the slurm header that I left in case one wants to run it there instead. But run times are so low that it's not worth going into the queue.

# $1 is the fastq file
# $2 is the fasta file
# $3 is the batch identifier

# ml build-env/2020 samtools/1.10-foss-2018b bwa/0.7.17-patch-1-foss-2018b
name=`basename $1 .fq.gz`

# bwa parameters
# mem = run bwa mem (maximum exact match algorithm)
# -t = number of threads [1] 
# -k = Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19] 
# -a = output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments. 
mkdir -p indexes bwa_logs bam idxstats
bwa mem -t 8 -k 15 -a indexes/${2} $1 2>bwa_logs/${name}.log | samtools view -b -o bam/${name}_${3}.bam -
samtools sort -@ 8 -o bam/${name}_${3}_sorted.bam bam/${name}_${3}.bam
samtools index bam/${name}_${3}_sorted.bam
samtools idxstats bam/${name}_${3}_sorted.bam > idxstats/${name}_${3}_sorted.bam.idxstats
rm bam/${name}_${3}.bam
