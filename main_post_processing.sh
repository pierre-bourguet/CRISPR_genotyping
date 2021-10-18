#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-00:10:00
#SBATCH --qos=short
#SBATCH --output=logs/%A_%a_%x.out
#SBATCH --error=logs/%A_%a_%x.err

id=$1

echo -e "$(date) .. Creating indel frequency tables, high coverage mutations..." # compile the data into one large table: selecting the mutation with the highest coverage, and extract frequency. The reason is that high coverage / frequency deletions are often followed by low coverage & even-higher-frequency mutations
c=0 && for a in mutations/*mutations.tsv ; do
    if [ $c = 0 ]; then echo > header && cut -f1 $a | uniq > tmp && echo "sum" >> tmp && tr "\n" "\t" < tmp > ${id}_frequency_of_top_coverage_indels.tsv && echo "" >> ${id}_frequency_of_top_coverage_indels.tsv && c=1 ; fi # row names
    touch ${a}_tmp
    for gene in `cut -f1 $a | sort | uniq` ; do
        grep $gene $a | sort -nrk4,4 -nrk6,6 | head -1 >> ${a}_tmp # mutation with the highest coverage (not frequency!) is reported. In case of equal coverage, the highest frequency is selected
    done    
	cut -f6 ${a}_tmp | tr "\n" "\t" | cat ${id}_frequency_of_top_coverage_indels.tsv - > tmp && echo `cut -f6 ${a}_tmp | sed '/^$/d' | paste -sd+ | bc` >> tmp && mv tmp ${id}_frequency_of_top_coverage_indels.tsv && rm ${a}_tmp
    echo -e "${a%_mutations.tsv}" | sed 's/mutations\///g' >> header
done
nbcol=$(echo `head -1 ${id}_frequency_of_top_coverage_indels.tsv | wc | cut -f13 -d " "` "+1" | bc) # put the number of columns in a variable
paste <(tail -n+2 header) <(tail -n+2 ${id}_frequency_of_top_coverage_indels.tsv) | sort -nrk${nbcol},${nbcol} -t$'\t' | cat <(head -1 ${id}_frequency_of_top_coverage_indels.tsv | paste <(head -1 header) -) - > tmp && mv tmp ${id}_frequency_of_top_coverage_indels.tsv
# same compiling but this time total frequency of mutations that passed the treshold
echo -e "$(date) .. Creating indel frequency tables, total mutation frequency..."
c=0 && for a in mutations/*mutations.tsv ; do
    if [ $c = 0 ]; then echo > header && cut -f1 $a | sort | uniq > tmp && echo "sum" >> tmp && tr "\n" "\t" < tmp > ${id}_total_frequency_of_indels.tsv && echo "" >> ${id}_total_frequency_of_indels.tsv && c=1 ; fi # row names
    touch ${a}_tmp
    for gene in `cut -f1 $a | sort | uniq` ; do
        grep $gene $a | awk -v tot=`grep $gene $a | cut -f6 | paste -sd+ | bc` 'BEGIN{FS=OFS="\t"};{print $1,$2,$3,$4,$5,tot}' | head -1 >> ${a}_tmp # mutation with the highest coverage (not frequency!) is reported
    done    
	cut -f6 ${a}_tmp | tr "\n" "\t" | cat ${id}_total_frequency_of_indels.tsv - > tmp && echo `cut -f6 ${a}_tmp | sed '/^$/d' | paste -sd+ | bc` >> tmp && mv tmp ${id}_total_frequency_of_indels.tsv && rm ${a}_tmp
    echo -e "${a%_mutations.tsv}" | sed 's/mutations\///g' >> header
done
nbcol=$(echo `head -1 ${id}_total_frequency_of_indels.tsv | wc | cut -f13 -d " "` "+1" | bc) # put the number of columns in a variable
paste <(tail -n+2 header) <(tail -n+2 ${id}_total_frequency_of_indels.tsv) | sort -nrk${nbcol},${nbcol} -t$'\t' | cat <(head -1 ${id}_total_frequency_of_indels.tsv | paste <(head -1 header) -) - > tmp && mv tmp ${id}_total_frequency_of_indels.tsv

echo -e "$(date) .. Cleaning up..."
mkdir -p $id && rm -r idxstats header
mv bam indexes bwa_logs mpileup mutations *tsv $id
