#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 56  # cores (note TACC allocates whole nodes)
#SBATCH -t 0-09:00  # runtime in D-HH:MM
#SBATCH -p small  # partition to submit to
##SBATCH --mem=192G  # total mem (note TACC allocates whole nodes; using this argument causes job to fail)
#SBATCH -J megahit_60m
#SBATCH -o slurm/megahit_60m.out
#SBATCH -e slurm/megahit_60m.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=
#SBATCH --dependency=afterok:4219445 # job can begin after specified jobs have completed with exit code zero

source activate snakemake

snakemake -j 56 --use-conda --rerun-incomplete --unlock out/megahit/60m/final.contigs.fa
snakemake -j 56 --use-conda --rerun-incomplete out/megahit/60m/final.contigs.fa
