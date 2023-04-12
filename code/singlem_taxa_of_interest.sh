#!/bin/bash

for file in out/singlem/*.tsv
do
    echo $file
    grep 'Stenotrophomonas' $file
    grep 'Parvibaculum' $file
    echo
done > out/singlem/taxa_of_interest.txt

# spare code
# grep 'Stenotrophomonas' all_good_samples.otu_table.tsv | awk '{print $1,$2,$4}' >singlem_results_stenotrophomonas.txt

for file in out/singlem/*.tsv
do
    grep 'Xanthomonadaceae$' $file
    echo
done > out/singlem/taxa_of_interest_Xanthomonadaceae.txt

for file in out/singlem/*.tsv
do
    grep 'Parvibaculaceae$' $file
    echo
done > out/singlem/taxa_of_interest_Parvibaculaceae.txt


