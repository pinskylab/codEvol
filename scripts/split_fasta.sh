#!/bin/bash

# inspired by https://www.biostars.org/p/105388/

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fasta
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < gadMor2.fasta
