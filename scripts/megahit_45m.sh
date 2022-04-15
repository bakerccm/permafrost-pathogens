#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 56  # cores (note TACC allocates whole nodes)
#SBATCH -t 0-09:00  # runtime in D-HH:MM
#SBATCH -p small  # partition to submit to
##SBATCH --mem=192G  # total mem (note TACC allocates whole nodes; using this argument causes job to fail)
#SBATCH -J megahit_45m
#SBATCH -o slurm/megahit_45m.out
#SBATCH -e slurm/megahit_45m.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

source activate snakemake

snakemake -j 56 --use-conda --rerun-incomplete --unlock out/megahit/45m/final.contigs.fa
snakemake -j 56 --use-conda --rerun-incomplete out/megahit/45m/final.contigs.fa
