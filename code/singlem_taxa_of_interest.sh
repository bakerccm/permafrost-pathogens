#!/bin/bash

for file in out/singlem/*.tsv
do
    echo $file
    grep 'Stenotrophomonas' $file
    grep 'Parvibaculum' $file
done > out/singlem/taxa_of_interest.txt

# spare code
# grep 'Stenotrophomonas' all_good_samples.otu_table.tsv | awk '{print $1,$2,$4}' >singlem_results_stenotrophomonas.txt

