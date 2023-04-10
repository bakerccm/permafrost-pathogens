#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 64  # cores - many cores have 24 cores on Premise but some have more
#SBATCH -t 3-00:00  # runtime in D-HH:MM
#SBATCH -p shared  # partition to submit to
#SBATCH --mem=950G  # total mem
#SBATCH --nodelist=node117 # comma separated list of node names to use`
#SBATCH -J snakemake
#SBATCH -o slurm/snakemake-%j.out
#SBATCH -e slurm/snakemake-%j.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=

##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

eval "$(conda shell.bash hook)"
conda activate snakemake-7.25.0

snakemake -j 64 --use-conda --rerun-incomplete --unlock out/megahit/{35m,45m,60m,83m,NT}/final.contigs.fa
snakemake -j 64 --use-conda --rerun-incomplete out/megahit/{35m,45m,60m,83m,NT}/final.contigs.fa

# snakemake -j 56 --use-conda out/megahit/35m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/45m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/60m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/83m/final.contigs.fa
# snakemake -j 56 --use-conda out/megahit/NT/final.contigs.fa

# snakemake -j 8 --use-conda some_rule --forcerun upstream_rule

# snakemake -j 1 --use-conda -rerun-incomplete --unlock some_rule
# snakemake -j 8 --use-conda -rerun-incomplete some_rule

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

