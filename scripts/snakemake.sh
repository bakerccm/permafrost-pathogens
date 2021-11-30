#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 56  # cores
#SBATCH -t 0-02:00  # runtime in D-HH:MM
#SBATCH -p development  # partition to submit to
#SBATCH --mem=192G  # memory pool for all cores (see also --mem-per-cpu)
#SBATCH -J snakemake
#SBATCH -o slurm/snakemake.out
#SBATCH -e slurm/snakemake.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

conda activate snakemake

snakemake -j 56 --use-conda cutadapt_multiQC
