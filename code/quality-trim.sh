#!/bin/bash
#$ -q TELOMER,UI-HM
#$ -pe smp 32
#$ -cwd
#$ -r y
# -j y
# -o /dev/null
#$ -l datacenter=LC


module unload jdk
module load openjdk/11.0.2

dir="./"


for i in ${dir}*.fastq.gz; do
    r1=$i
    r2=${i/R1/R2}
    o1=${r1%.fastq.gz}.p.fq.gz
    o1_up=${r1%.fastq.gz}.up.fq.gz
    o2=${r2%.fastq.gz}.p.fq.gz
    o2_up=${r2%.fastq.gz}.up.fq.gz

    java -jar ~/bin/trimmomatic-0.39.jar PE -threads 32 -phred33 ${r1} ${r2} ${o1} ${o1_up} ${o2} ${o2_up} ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 $


done


for i in ${dir}*.p.fq.gz; do
    o3=${i%.fq.gz}.fin.fastq
    prinseq++ -threads 32 -fastq ${i} -lc_dust=0.07 -out_good ${o3} -out_bad /dev/null

done

