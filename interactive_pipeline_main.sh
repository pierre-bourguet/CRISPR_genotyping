#!/usr/bin/env bash

# arguments
indel=5 # size of indels that will be tested
threshold=0.1 # frequency threshold (from 0 to 1) to report mutations

# import arguments
while getopts ":q:f:n:" options; do
    case "${options}" in
        q) fq_dir=${OPTARG%/};;
        f) fasta=${OPTARG};;
        n) id=${OPTARG};;
    esac
done

if [ "$fq_dir" = "" ]; then
    echo -e "-q is empty"
elif [ "$fasta" = "" ]; then
    echo -e " -f is empty"
elif [ "$id" = "" ]; then
    echo -e " -n is empty"
fi

echo -e "$(date) .. Running pipeline for $id using fastq files in $fq_dir directory against $fasta reference..."

ml build-env/2020 samtools/1.10-foss-2018b bwa/0.7.17-patch-1-foss-2018b
fasta_basename=`basename $fasta`

echo -e "$(date) .. Constructing index..."
bwa index $fasta
mkdir -p indexes && mv ${fasta}?* indexes && cp ${fasta} indexes/

echo -e "$(date) .. Mapping with bwa..."
for i in ${fq_dir}/*gz ; do bash interactive_pipeline_run_bwa_mem.sh $i $fasta_basename $id ; done

echo -e "$(date) .. Create idxstat summary file..."
cd idxstats && c=0 && for i in * ; do
    if [ $c = 0 ]; then cut -f1,2 $i > ../${id}_idxstats.tsv && touch ../header && c=1 ; fi 
    name=`echo $i | sed "s/${id}_sorted.bam.idxstats// ; s/${id}_multiplex_CRISPR_//"` && cut -f3 $i | paste ../${id}_idxstats.tsv - > tmp && mv tmp ../${id}_idxstats.tsv && echo -e $name >> ../header ; done
cd .. && tr "\n" "\t" < header | sed 's/^/\t\t/ ; s/$/\n/' | cat - ${id}_idxstats.tsv > tmp && mv tmp ${id}_idxstats.tsv && rm header

echo -e "$(date) .. Running samtools mpileup..."
# create pileup file to summarize mismatches
mkdir -p mpileup mutations
for i in bam/*bam ; do
    output=`echo "${i##*/}" | sed "s/${id}_multiplex_CRISPR_// ; s/_${id}_sorted.bam//"`
    samtools mpileup -f indexes/${fasta_basename} -o mpileup/${output}_mpileup.tsv $i 
done

echo -e "$(date) .. Counting frequency of indels in mpileup file..."
for a in mpileup/*mpileup.tsv ; do
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
    name=`basename $a .tsv`
    sort -k1,1 -k2,2n ${a}.tmp > tmp && mv tmp mutations/${name%_mpileup}_mutations.tsv && rm ${a}.tmp
done

echo -e "$(date) .. Creating indel frequency tables..."
# compile the data into one large table: selecting the mutation with the highest coverage, and extract frequency. The reason is that high coverage / frequency deletions are often followed by low coverage & even-higher-frequency mutations

c=0 && for a in mutations/*mutations.tsv ; do
    if [ $c = 0 ]; then echo > header && cut -f1 $a | sort | uniq > tmp && echo "sum" >> tmp && tr "\n" "\t" < tmp > ${id}_frequency_of_top_coverage_indels.tsv && echo "" >> ${id}_frequency_of_top_coverage_indels.tsv && c=1 ; fi # row names
    touch ${a}_tmp
    for gene in `cut -f1 $a | sort | uniq` ; do
        grep $gene $a | sort -nrk4,4 | head -1 >> ${a}_tmp # mutation with the highest coverage (not frequency!) is reported
    done    
	cut -f6 ${a}_tmp | tr "\n" "\t" | cat ${id}_frequency_of_top_coverage_indels.tsv - > tmp && echo `cut -f6 ${a}_tmp | sed '/^$/d' | paste -sd+ | bc` >> tmp && mv tmp ${id}_frequency_of_top_coverage_indels.tsv && rm ${a}_tmp
    echo -e "${a%_mutations.tsv}" | sed 's/mutations\///g' >> header
done
echo " `wc -l header | cut -f1 -d " " `" "+ 1 " | bc > nb_col && nbcol=`cat nb_col` # put the number of columns in a variable
sort -rk${nbcol},${nbcol} ${id}_frequency_of_top_coverage_indels.tsv | paste header - > tmp && mv tmp ${id}_frequency_of_top_coverage_indels.tsv

# same compiling but this time total frequency of mutations
c=0 && for a in mutations/*mutations.tsv ; do
    if [ $c = 0 ]; then echo > header && cut -f1 $a | sort | uniq > tmp && echo "sum" >> tmp && tr "\n" "\t" < tmp > ${id}_total_frequency_of_indels.tsv && echo "" >> ${id}_total_frequency_of_indels.tsv && c=1 ; fi # row names
    touch ${a}_tmp
    for gene in `cut -f1 $a | sort | uniq` ; do
        grep $gene $a | awk -v tot=`grep $gene $a | cut -f6 | paste -sd+ | bc` 'BEGIN{FS=OFS="\t"};{print $1,$2,$3,$4,$5,tot}' | head -1 >> ${a}_tmp # mutation with the highest coverage (not frequency!) is reported
    done    
	cut -f6 ${a}_tmp | tr "\n" "\t" | cat ${id}_total_frequency_of_indels.tsv - > tmp && echo `cut -f6 ${a}_tmp | sed '/^$/d' | paste -sd+ | bc` >> tmp && mv tmp ${id}_total_frequency_of_indels.tsv && rm ${a}_tmp
    echo -e "${a%_mutations.tsv}" | sed 's/mutations\///g' >> header
done
sort -rk${nbcol},${nbcol} ${id}_total_frequency_of_indels.tsv | paste header - > tmp && mv tmp ${id}_total_frequency_of_indels.tsv

echo -e "$(date) .. Cleaning up..."
mkdir -p $id && rm -r idxstats header nb_col
mv bam indexes bwa_logs mpileup mutations *tsv $id