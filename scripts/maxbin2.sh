#!/bin/bash
#SBATCH -N 2  # nodes
#SBATCH -n 48  # cores (note TACC allocates whole nodes)
#SBATCH -t 2-00:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
##SBATCH --mem=192G  # total memory
#SBATCH -J maxbin2
#SBATCH -o slurm/maxbin2-%j.out
#SBATCH -e slurm/maxbin2-%j.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

eval "$(conda shell.bash hook)"
conda activate snakemake-7.25.0

snakemake -j 48 --use-conda --rerun-incomplete --unlock all_maxbin2

snakemake -j 48 --use-conda --rerun-incomplete all_maxbin2

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

