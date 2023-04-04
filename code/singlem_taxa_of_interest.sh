#!/bin/bash

for file in out/singlem/*.tsv
do
    echo $file
    grep 'Stenotrophomonas' $file
    grep 'Parvibaculum' $file
done > out/singlem/taxa_of_interest.txt

