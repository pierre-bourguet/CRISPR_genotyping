#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-00:05:00
#SBATCH --qos=short
#SBATCH --output=logs/%A_%a_%x.out
#SBATCH --error=logs/%A_%a_%x.err

# $1 is the fastq directory
# $2 is the fasta file
# $3 is the batch identifier
# $4 is the threshold frequency for reporting mutations
# $5 is the maximum tested indel size

# ml build-env/2020 samtools/1.10-foss-2018b bwa/0.7.17-patch-1-foss-2018b
fq_file=`ls $1 | sed -n ${SLURM_ARRAY_TASK_ID}p`
# name=`basename $1 .fq.gz`
name=`basename $fq_file .fq.gz`
id=$3
threshold=$4 # frequency threshold (from 0 to 1) to report mutations
indel=$5 # size of indels that will be tested

echo -e "$(date) .. Running pipeline for $3 using fastq files in $1 directory with file $fq_file against $2 reference..."

#### bwa mapping
# bwa parameters
# mem = run bwa mem (maximum exact match algorithm)
# -t = number of threads [1] 
# -k = Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19] 
# -a = output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
echo -e "$(date) .. Mapping with bwa..."
mkdir -p bwa_logs bam idxstats mpileup mutations
bwa mem -t 2 -k 15 -a indexes/${2} ${1}/${fq_file} 2>bwa_logs/${name}.log | samtools view -b -o bam/${name}_${id}.bam -
echo -e "$(date) .. Sorting & indexing bam files..."
samtools sort -@ 2 -o bam/${name}_${id}_sorted.bam bam/${name}_${id}.bam
samtools index bam/${name}_${id}_sorted.bam
echo -e "$(date) .. Running idxstats..."
samtools idxstats bam/${name}_${id}_sorted.bam > idxstats/${name}_${id}_sorted.bam.idxstats
rm bam/${name}_${id}.bam

# create pileup file to summarize mismatches
echo -e "$(date) .. Running samtools mpileup..."
bam=bam/${name}_${id}_sorted.bam
output=`echo ${bam##*/} | sed 's/_sorted.bam//'`
samtools mpileup -f indexes/${2} -o mpileup/${output}_mpileup.tsv $bam 

echo -e "$(date) .. Counting frequency of indels in mpileup file & reporting mutations above ${threshold} frequency..."
a=mpileup/${output}_mpileup.tsv
rm -f ${a}.tmp && touch ${a}.tmp
for gene in `cut -f1 $a | sort | uniq` ; do
    for i in A T G C ; do 
        for j in \- \+ ; do
            for k in $(seq ${indel}) ; do
                grep $gene $a | awk -F'\t' -v nt="${j}${k}${i}" '{print $0, "\t" gsub(nt,"")}' | awk -F'\t' 'BEGIN{OFS="\t"} $4 > 0 {$8= $7 / $4}1' | awk -v nt="${j}${k}${i}" -v t="$threshold" 'BEGIN{OFS="\t"} $8>t {print $1,$2,$3,$4,nt,$8}' >> ${a}.tmp
            done
        done
    done
if [ `grep $gene ${a}.tmp | wc -l | cut -f1 | cut -f1 -d " "` = 0 ]; then echo -e "${gene} `printf '\t%.0s' {1..5}`" >> ${a}.tmp ; fi # echo an empty line if no mutations above the threshold was found
done
name=`basename $a .tsv | sed 's/_sorted.bam//'`
sort -k1,1 -k2,2n ${a}.tmp > tmp && mv tmp mutations/${name%_mpileup}_mutations.tsv && rm ${a}.tmp

