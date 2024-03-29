#!/bin/bash
#SBATCH -N 5 # nodes
##SBATCH -n 120  # cores
#SBATCH -t 0-03:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
##SBATCH --mem=192G  # total memory
#SBATCH -J all_prokka
#SBATCH -o slurm/all_prokka-%j.out
#SBATCH -e slurm/all_prokka-%j.err
#SBATCH --nodelist=node136,node137,node138,node139,node140 # comma separated list of node names to use
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=bakerccm@gmail.com

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

eval "$(conda shell.bash hook)"
conda activate snakemake-7.25.0

snakemake -j 1 --use-conda --rerun-incomplete --unlock all_prokka_donefiles
snakemake -j 120 --use-conda --rerun-incomplete all_prokka_donefiles

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

