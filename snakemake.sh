#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 8    # cores
#SBATCH --mem=64G    # memory per node
##SBATCH --mem-per-cpu=6G    # memory per cpu
#SBATCH -p priority    # partition
#SBATCH -t 1-00:00    # runtime d-hh:mm
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err
#SBATCH --mail-user=xxxxxxxx@xxxxx.com
#SBATCH --mail-type=END
##SBATCH --dependency=afterOK:9999999

# module load snakemake/3.12.0

# use snakemake environment with snakemake 5.10.0 installed (env modules on cluster only have snakemake/3.12.0)
# note use of source activate (conda activate became preferred from conda v4.4)
module load conda2/4.2.13
source activate snakemake

snakemake --use-conda -j 8 some_rule

# snakemake -j 8 --use-conda some_rule --forcerun upstream_rule

# snakemake -j 1 --use-conda -rerun-incomplete --unlock some_rule
# snakemake -j 8 --use-conda -rerun-incomplete some_rule
