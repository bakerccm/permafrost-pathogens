# Shotgun metagenomic data for permafrost pathogens project

The data for this project are not included in the code repository and should be downloaded to this folder in order to reproduce the pipeline results.

The data for this project were originally published with [Barbato **et al.** (2022)](#citations).

The data are publicly available and may be downloaded from the NCBI.

## Files in this folder

`fastq` - raw data should be saved into this folder (unless a different location is specified in the [config file](/config/config.yaml)

md5sums.txt - list of md5 hashes for the uncompressed fastq files; note that these were generated with 6 threads in fasterq-dump (different thread counts will affect the md5 hash); it is probably more useful and reliable to use vdb-validate to check the integrity of the data download

readme.md - is this file.

download.sh - provides download commands. you may need to modify these commands, depending on how you are downloading the data


## Citations

Barbato, RA, RM Jones, TA Douglas, SJ Doherty, K Messan, KL Foley, EJ Perkins, AK Thurston and N Garcia-Reyero. Not all permafrost microbiomes are created equal: Influence of permafrost thaw on the soil microbiome in a laboratory incubation study (2022). *Soil Biology and Biochemistry* 167:108605 [doi:10.1016/j.soilbio.2022.108605](https://doi.org/10.1016/j.soilbio.2022.108605) [[https://www.sciencedirect.com/science/article/pii/S0038071722000621]](https://www.sciencedirect.com/science/article/pii/S0038071722000621)
