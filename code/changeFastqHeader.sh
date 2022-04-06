#!/bin/bash

if [ $# != "2" ] || [ $1 = "-h" ] || [ $1 = "--help" ]; then
        echo "$0 fqin fqout"
        echo "Rewrite fastq header to use output target name"
        echo "fqin: fasta input file"
        echo "fqout: fasta output file"
        echo ""
        exit 1
fi

fqin=$1
fqout=$2

# Grab just the filename, minus extension
file=$(basename $fqout)
filename=${file%_R*}

awkscript=\
'{
        if (NR % 4 == 1) {
                print "@" name "_" ++i
        } else if (NR % 4 == 3) {
                print "+" name "_" ++j
        } else {
                print $0
        }
}'

awk -v "name=$filename" "$awkscript" $fqin > $fqout
