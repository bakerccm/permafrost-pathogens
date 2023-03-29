#!/bin/bash
#SBATCH -N 1             # request one whole node (should have 24 cores)
#SBATCH -n 6             # request eight cores
#SBATCH --mem=16G        # request 16GB memory on the node
#SBATCH -t 0-01:00       # runtime in D-HH:MM
#SBATCH -p shared        # partition to submit to
#SBATCH -J fasterq-dump  # job name
#SBATCH -o slurm/fasterq-%j.out
#SBATCH -e slurm/fasterq-%j.err
##SBATCH --mail-type=BEGIN,END    # notifications: BEGIN,END,FAIL,ALL

# get accession number from command line

accession=$1

source activate SRAtools

cd ~/permafrost-pathogens/data/fastq

# for accession in `tail -n +40 ../../metadata/SraAccList.txt | head -n 1`
# do

    echo ${accession}
    # get fastq files; 6 threads is the default
        fasterq-dump --threads 6 --temp . ${accession}/${accession}.sra
    # compress fastq files
        gzip ${accession}.sra_1.fastq &
        gzip ${accession}.sra_2.fastq &
        wait
    # remove .sra files
    #   rm ${accession}/${accession}.sra && rmdir ${accession}
   echo ${accession} complete

# done

sleep 5
sacct -j $SLURM_JOBID --format=JobID,JobName%24,CPUTime,Elapsed,MaxRSS --units=G

