#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 24  # cores - many cores have 24 cores on Premise but some have more
#SBATCH -t 1-00:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
#SBATCH --mem=500G  # total mem
##SBATCH --nodelist=node117 # comma separated list of node names to use
#SBATCH -J singlem
#SBATCH -o slurm/singlem-%j.out
#SBATCH -e slurm/singlem-%j.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

# use snakemake conda environment with snakemake 6.4.1 installed
# note use of conda activate (preferred over source activate from conda v4.4 onwards)
# conda activate snakemake

source activate snakemake-6.4.1

snakemake -j 24 --use-conda --rerun-incomplete --unlock all_singlem
snakemake -j 24 --use-conda --rerun-incomplete all_singlem

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

