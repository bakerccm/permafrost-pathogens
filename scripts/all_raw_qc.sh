#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 58  # cores
#SBATCH -t 0-01:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
##SBATCH --mem=192G  # total mem
#SBATCH -J all_raw_qc
#SBATCH -o slurm/all_raw_qc.out
#SBATCH -e slurm/all_raw_qc.err
#SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bakerccm@gmail.com

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

# use snakemake conda environment with snakemake 6.4.1 installed
# note use of conda activate (preferred over source activate from conda v4.4 onwards)
# conda activate snakemake

source activate snakemake

snakemake -j 58 --use-conda all_raw_qc
