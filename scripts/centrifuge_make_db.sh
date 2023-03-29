#!/bin/bash
#SBATCH -N 1  # nodes
#SBATCH -n 56  # cores (note TACC allocates whole nodes)
#SBATCH -t 2-00:00  # runtime in D-HH:MM
#SBATCH -p small  # partition to submit to
##SBATCH --mem=192G  # total mem (note TACC allocates whole nodes; using this argument causes job to fail)
#SBATCH -J centrifuge_make_db
#SBATCH -o slurm/centrifuge_make_db-%j.out
#SBATCH -e slurm/centrifuge_make_db-%j.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=
##SBATCH --dependency=after:jobid[:jobid...] job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] job can begin after specified jobs have failed

source activate centrifuge

# cd $WORK/permafrost-pathogens/out/anvio/centrifuge

cd $SCRATCH/centrifuge

make THREADS=52 DONT_DUSTMASK=1 p_compressed

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

