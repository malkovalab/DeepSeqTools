#!/bin/bash
#$ -q TELOMER,UI-HM
#$ -pe smp 56
#$ -cwd
#$ -r y
# -j y
# -o /dev/null
#$ -l datacenter=LC


dir="raw_reads/"
ref="ref/lys2InsH.fa"

#for batch processing


for i in ${dir}*R1*.p.fin.fastq; do
    r1=$i
    r1_NH=${r1}.NH.fq
    ./changeFastqHeader.sh $r1 $r1_NH

    r2=${i/R1/R2}
    r2_NH=${r2}.NH.fq
    ./changeFastqHeader.sh $r2 $r2_NH

    OUT=${i%_R1_001.fastq.gz}
    OUT=${OUT#$dir}
    echo $OUT
    bwa mem -M -t 56 $ref $r1_NH $r2_NH > $OUT.PP.sam
    samtools sort -@56 -m 2G -O bam -T temp.sort -o $OUT.PP.sorted.bam $OUT.PP.sam

done
