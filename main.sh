#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-00:02:00
#SBATCH --qos=short
#SBATCH --output=logs/%A_%a_%x.out
#SBATCH --error=logs/%A_%a_%x.err

# arguments
threshold=0.1 # frequency threshold (from 0 to 1) to report mutations
indel=30 # size of indels that will be tested

# import arguments
while getopts ":q:f:n:" options; do
    case "${options}" in
        q) fq_dir=${OPTARG%/};;
        f) fasta=${OPTARG};;
        n) id=${OPTARG};;
    esac
done

echo -e "$(date) .. Running pipeline for $id using fastq files in $fq_dir directory against $fasta reference..."

ml build-env/2020 samtools/1.10-foss-2018b bwa/0.7.17-patch-1-foss-2018b
fasta_basename=`basename $fasta`

echo -e "$(date) .. Constructing index..."
bwa index $fasta
mkdir -p indexes && mv ${fasta}?* indexes && cp ${fasta} indexes/

echo -e "$(date) .. Submitting job arrays to map with bwa & count indels for each sample..."
nb_sample=` ls ${fq_dir}/*gz | wc -l `
jid1=$(sbatch --parsable --array=1-${nb_sample} run_bwa_mem.sh $fq_dir $fasta_basename $id $threshold $indel)
jid2=$(sbatch --parsable --depend=afterany:$jid1 main_post_processing.sh $id)
