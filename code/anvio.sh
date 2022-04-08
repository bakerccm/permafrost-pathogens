################################
# install anvio

# installation involves manual step
# note sure if it can be done through conda in snakemake (i.e. without manual intervention), but if not then how?
# see https://stackoverflow.com/questions/59107413/activating-existing-conda-enviornments-in-snakemake

#following anvio install instructions from https://anvio.org/install/ except with use of yaml file for first bit
mamba env create -f envs/anvio.yaml
# just kidding, it fails with error
#      ModuleNotFoundError: No module named '_sysconfigdata_x86_64_conda_linux_gnu'
# ah, see the fix for that at https://anvio.org/install/ :
# see also discussion here https://github.com/merenlab/anvio/issues/1839
 cd ~/miniconda3/envs/anvio/lib/python3.6
 mv _sysconfigdata_x86_64_conda_cos6_linux_gnu.py _sysconfigdata_x86_64_conda_linux_gnu.py
 cd $pathogens # or whatever
# can this pip part also be done in the yaml file?
curl -L https://github.com/merenlab/anvio/releases/download/v7.1/anvio-7.1.tar.gz \
        --output anvio-7.1.tar.gz
conda activate anvio
# ok that seems to have worked

################################
# mapping to read assembly

# see https://astrobiomike.github.io/metagenomics/metagen_anvio

# reformat deflines in FASTA file containing contigs from megahit
conda activate anvio
anvi-script-reformat-fasta "out/megahit/35m-t0-R1/final.contigs.fa" -o "out/megahit/35m-t0-R1/contigs.fa" -l 1000 --simplify-names

# the remainder of this part can be snakemakified with conda as usual

# create an index of co-assembly
# for each co-assembly:
conda activate bowtie2
mkdir -p "out/assemblies/35m-t0-R1" # at level of assembly or co-assembly
bowtie2-build "out/megahit/35m-t0-R1/contigs.fa" "out/assemblies/35m-t0-R1/assembly" # note contigs.fa is the version with the cleaned deflines

# run mapping
# took about 20 min to do this on idev for this one sample
# for each sample, and for each co-assembly you want to map it to:
mkdir -p "out/mapping/35m-t0-R1" # folder reflects co-assembly being mapped to
bowtie2 -x "out/assemblies/35m-t0-R1/assembly" -q \
-1 "out/bbduk_noPhiX_fastuniq/35m-t0-R1_R1.fastq.gz" \
-2 "out/bbduk_noPhiX_fastuniq/35m-t0-R1_R2.fastq.gz" \
--no-unal -p 8 -S "out/mapping/35m-t0-R1/35m-t0-R1.sam" # folder name reflects co-assembly being mapped to; file name reflects sample being mapped

# convert sam to bam
conda activate samtools # (actually, no need to do this , samtools is already in the anvio conda environment)
# 5 min with one sample
samtools view -b -o "out/mapping/35m-t0-R1/35m-t0-R1_raw.bam" "out/mapping/35m-t0-R1/35m-t0-R1.sam"

# sort and index bam file
# a couple of minutes on idev
samtools sort -o "out/mapping/35m-t0-R1/35m-t0-R1.bam" "out/mapping/35m-t0-R1/35m-t0-R1_raw.bam"
samtools index "out/mapping/35m-t0-R1/35m-t0-R1.bam"

################################
# anvio part
# see https://astrobiomike.github.io/metagenomics/metagen_anvio#mapping-our-reads-to-the-assembly-they-built

# generate anvi’o contigs database from co-assembly fasta
mkdir -p out/anvio
conda activate anvio # assuming we were in the samtools environment from above
anvi-gen-contigs-database -f "out/megahit/35m-t0-R1/contigs.fa" -o "out/anvio/contigs.db" -n "metagenome-35m-t0-R1"

# scan for commonly used bacterial single-copy genes to help estimate
# genome completeness/redundancy in real-time as we bin our contigs later
anvi-run-hmms -c out/anvio/contigs.db -T 4 # tutorial adds -I Campbell_et_al but this profile is not available in my anvio install - presumably need to get from the paper

# use NCBI COGs for functional annotation
# COGs = clusters of orthologous genes
# homologs = genes with shared ancestry
# of which:
#     orthologs = genes that diverge due to genome duplication (i.e. speciation)
#     paralogs = genes that diverge due to gene duplication
#     xenologs = genes that share ancestry due to HGT
anvi-setup-ncbi-cogs --cog-data-dir out/anvio/COGs -T 4 # --just-do-it
anvi-run-ncbi-cogs --cog-data-dir out/anvio/COGs -c out/anvio/contigs.db  -T 4

# assign taxonomy with centrifuge
# prepare for centrifuge
anvi-get-sequences-for-gene-calls -c "out/anvio/contigs.db" -o "out/anvio/gene_calls.fa"
# install centrifuge and activate environment (actually, no need, centrifuge is already in the anvio conda environment)
mamba env create -f envs/centrifuge.yaml # installs centrifuge, jellyfish and mummer
conda activate centrifuge
# prepare taxonomy database
# follow links from https://github.com/infphilo/centrifuge on pre-built indices
# --> http://www.ccb.jhu.edu/software/centrifuge/
mkdir -p out/anvio/centrifuge
cd out/anvio/centrifuge
wget https://raw.githubusercontent.com/infphilo/centrifuge/master/indices/Makefile
# currently 415 archaeal, 25823 bacterial and xxxx viral genomes
# make THREADS=48 p+h+v    # bacterial, human, and viral genomes [~12G]
make THREADS=52 p_compressed    # bacterial genomes compressed at the species level [~4.2G]
make THREADS=52 DONT_DUSTMASK=1 p_compressed # use DONT_DUSTMASK because of missing NBCI tools
# make THREADS=52 p_compressed+h+v    # combination of the two above [~8G]
cd ../../..

# run centrifuge
centrifuge -f -x out/anvio/centrifuge/p+h+v "out/anvio/gene_calls.fa" -S "out/anvio/centrifuge_hits.tsv" -p 16

# import centrifuge results into anvio contigs database
anvi-import-taxonomy-for-genes -c out/anvio/contigs.db -i out/anvio/centrifuge_report.tsv out/anvio/centrifuge_hits.tsv -p centrifuge

# profile samples
# output is in out/mapping/35m-t0-R1/s35m_t0_R1/PROFILE.db
anvi-profile -i "out/mapping/35m-t0-R1/35m-t0-R1.bam" -c "out/anvio/contigs.db" -T 8 --cluster-contigs # how to provide output location?
# I think you may not need to use --cluster-contigs if working with multiple samples?

# merge sample-level profiles when we are doing that
# anvi-merge */PROFILE.db -o merged_profile -c contigs.db
# anvi-merge "out/anvio/*/profile.db" -o "out/anvio/merged_profile" -c "out/anvio/contigs.db"

################################
# examine interactively

# ON FRONTERA

# download data for local examine, or see instructions here
# https://merenlab.org/2015/11/28/visualizing-from-a-server/
# for connecting using an SSH tunnel

# log on to Frontera with port forwarding
ssh -L 8080:localhost:8080 cbaker@frontera.tacc.utexas.edu
idev -t 1:30:00
cd $pathogens
conda activate anvio
anvi-interactive -p out/mapping/35m-t0-R1/s35m_t0_R1/PROFILE.db -c out/anvio/contigs.db --server-only -P 8080
# then open local brower and go to
http://localhost:8080

# LOCALLY

# move PROFILE.db and contigs.db to local folder:
cd /Users/chris/Documents/work/CRREL/permafrost-pathogens/anvio
conda activate anvio-7.1 # note slightly different environment name locally
anvi-interactive -p PROFILE.db -c contigs.db
