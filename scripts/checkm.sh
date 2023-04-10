#!/bin/bash
#SBATCH -N 1 # nodes
#SBATCH -n 24  # cores
#SBATCH -t 0-04:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
##SBATCH --mem=192G  # total memory
#SBATCH -J checkm
#SBATCH -o slurm/checkm-%j.out
#SBATCH -e slurm/checkm-%j.err
#SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=bakerccm@gmail.com

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

eval "$(conda shell.bash hook)"
conda activate snakemake-7.25.0

snakemake -j 1 --use-conda --rerun-incomplete --unlock out/maxbin2_checkm/{35m,45m,60m,83m,NT}.txt
snakemake -j 24 --use-conda --rerun-incomplete out/maxbin2_checkm/{35m,45m,60m,83m,NT}.txt

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

