#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 56  # cores
#SBATCH -t 0-03:00  # runtime in D-HH:MM
#SBATCH -p development  # partition to submit to
#SBATCH --mem-per-cpu=2G  # memory per CPU (see also --mem)
#SBATCH -J snakemake
#SBATCH -o slurm/snakemake.out
#SBATCH -e slurm/snakemake.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

# use snakemake conda environment with snakemake 6.4.1 installed
# note use of conda activate (preferred over source activate from conda v4.4 onwards)
conda activate snakemake

snakemake -j 56 --use-conda sickle_multiQC fastq_join_multiQC

# snakemake -j 8 --use-conda some_rule --forcerun upstream_rule

# snakemake -j 1 --use-conda -rerun-incomplete --unlock some_rule
# snakemake -j 8 --use-conda -rerun-incomplete some_rule
