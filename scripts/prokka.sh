#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 56  # cores (note TACC allocates whole nodes)
#SBATCH -t 0-01:00  # runtime in D-HH:MM
#SBATCH -p small  # partition to submit to
##SBATCH --mem=192G  # total mem (note TACC allocates whole nodes; using this argument causes job to fail)
#SBATCH -J prokka
#SBATCH -o slurm/prokka.out
#SBATCH -e slurm/prokka.err
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

snakemake -j 16 --use-conda out/maxbin2_prokka/35m/35m.{001..053}
snakemake -j 16 --use-conda out/maxbin2_prokka/45m/45m.{001..107}
snakemake -j 16 --use-conda out/maxbin2_prokka/60m/60m.{001..099}
snakemake -j 16 --use-conda out/maxbin2_prokka/83m/83m.{001..071}
snakemake -j 16 --use-conda out/maxbin2_prokka/NT/NT.{001..096}

