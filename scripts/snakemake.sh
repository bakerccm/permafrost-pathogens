#!/bin/bash
#SBATCH -N 2  # nodes
#SBATCH -n 48  # cores
#SBATCH -t 4-00:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
##SBATCH --mem=192G  # total mem
#SBATCH -J snakemake
#SBATCH -o slurm/snakemake.out
#SBATCH -e slurm/snakemake.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=
#SBATCH --gres=gpu:1

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

# use snakemake conda environment with snakemake 6.4.1 installed
# note use of conda activate (preferred over source activate from conda v4.4 onwards)
# conda activate snakemake

source activate snakemake

snakemake -j 48 --use-conda --rerun-incomplete --unlock out/megahit/{35m,83m}/final.contigs.fa
snakemake -j 48 --use-conda --rerun-incomplete out/megahit/{35m,83m}/final.contigs.fa

# snakemake -j 56 --use-conda out/megahit/35m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/45m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/60m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/83m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/NT/final.contigs.fa

# snakemake -j 8 --use-conda some_rule --forcerun upstream_rule

# snakemake -j 1 --use-conda -rerun-incomplete --unlock some_rule
# snakemake -j 8 --use-conda -rerun-incomplete some_rule
