#!/bin/bash
 
set -e
set -u
set -o pipefail
 
samples=${1:-./data/allbam.txt}
chromsize=${2:-./data/hg38.chrom.sizes}
size=${3:-3000}
COV_DIR="./data/genome_coverage"
CHR_DIR="./data/bam_chr_info"

mkdir -p ${COV_DIR}
mkdir -p ${CHR_DIR}

cat "$samples" | while read sample
do
    file_name=$(echo "$sample" | awk '{print $1}')
    file_ID=$(echo "$sample" | awk '{print $2}')
    echo ${file_name}
    echo ${file_ID}
    samtools view -H $sample | cut -f2 | grep SN | { while read top3 ; do echo ${top3:3} ; done } > ${CHR_DIR}/${file_ID}.txt
    bedtools makewindows -g $chromsize -w $size | bedtools sort -i -  -g $CHR_DIR/${file_ID}.txt | \
    bedtools intersect -b ${file_name} -a - -c -sorted -bed > ${COV_DIR}/${file_ID}.ReadsCoverage
done
