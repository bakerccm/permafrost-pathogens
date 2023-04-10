#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 24  # cores
#SBATCH -t 0-01:00  # runtime in D-HH:MM
#SBATCH -p shared # partition to submit to
##SBATCH --mem=500G  # total memory
#SBATCH -J metaquast
#SBATCH -o slurm/metaquast-%j.out
#SBATCH -e slurm/metaquast-%j.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

eval "$(conda shell.bash hook)"
conda activate snakemake-7.25.0

snakemake -j 24 --use-conda --rerun-incomplete out/metaquast/{35m,45m,60m,83m,NT}/report.txt

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

