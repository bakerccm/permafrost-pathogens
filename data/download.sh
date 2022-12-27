# Code for downloading data from NCBI to Premise

# Uses conda to install SRAtools to provide prefetch, vdb-validate and fasterq-dump.

# The environment module linuxbrew/colsa also provides prefetch and vdb-validate, but does
# not appear to provide fasterq-dump.

################################################################################################
# (1) Install conda

    # find latest version of miniconda at https://docs.conda.io/en/latest/miniconda.html and download
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh

    # run installer script (or could make it executable)
    bash Miniconda3-py39_4.9.2-Linux-x86_64.sh

    # add mamba installer to base environment
    conda install -c conda-forge mamba

################################################################################################
# (2) Install SRAtools from NCBI into a new conda environment

    # use mamba for the install since it's probably faster than conda
    mamba create -n SRAtools -c bioconda sra-tools


################################################################################################
# (3) Prepare list of accessions to download

    # See documentation on SRAtoolkit documentation here:
    #     https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc
    # Specifically on the prefetch command:
    #     https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=prefetch

    # Get SRA accession list file from NCBI site (SraAccList.txt) = one line per accession:

    #    SRR15048733
    #    SRR15048734
    #    SRR15048735
    #    ... etc

    # For this dataset, navigate to list of 58 SRA experiments at https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=542925
    # Then choose: Send to - File - Accession list
    # Downloaded file SraAccList.txt should have md5 hash dad41b5873fe3e85bf1e658f804ab54d

    # Transfer to cluster using scp or some other approach (e.g. cut and paste in nano)

    # Instructions below assume that the SRA accession list is saved at metadata/SraAccList.txt

################################################################################################
# (4) Prefetch data from NCBI

    # if helpful: start new screen session

        screen -S alaska

        screen -ls # check screen session; optional
        # type <ctrl>-a <ctrl>-d to detach
        screen -r alaska # to re-attach

    # download SRA files using prefetch from SRAtoolkit

        mkdir -p data/fastq # if folder doesn't already exist
        cd data/fastq

        conda activate SRAtools # activate environment to make SRAtools available

        # prefetch SRR15048733 # downloads a single accession
        # prefetch SRR15048735 SRR15048736 # downloads these two accessions
        prefetch $(<../../metadata/SraAccList.txt) # downloads all the accessions in the file by supplying the access numbers to the prefetch command as a space-separated list

        # Modify these commands as required

        # Note that this download proceeds serially through the accession list. It could easily take 12-18 hrs to complete,
        # so detach from screen while download is taking place. The screen session should stay active for the duration.

        # It is likely that the download could be done much faster in parallel. e.g. on the TACC, the entire
        # set of files could be downloaded in an interactive session in <30 min. On Premise this is
        # complicated by the use of a proxy server to provide internet connectivity on the compute nodes.

################################################################################################
# (5) Validate prefetched sra files

    # Get an interactive compute session
        salloc -p shared -t 0-02:00 -N 1 -n 1 --mem=16G
        # then connect to compute node using ssh

    conda activate SRAtools

    cd data/fastq
    tree
    # .
    # ├── SRR15048733
    # │   └── SRR15048733.sra
    # ├── SRR15048734
    # │   └── SRR15048734.sra
    # ├── SRR15048735
    # │   └── SRR15048735.sra
    # ...

    # vdb-validate is reasonably fast so a simple loop will do
    # but this could be parallelized if desired
    for accession in $(<../../metadata/SraAccList.txt)
    do
        echo ${accession}
        vdb-validate ${accession}/${accession}.sra
    done

    # SRR15048733
    # 2022-12-21T22:34:24 vdb-validate.2.10.0 info: Database 'SRR15048733.sra' metadata: md5 ok
    # 2022-12-21T22:34:24 vdb-validate.2.10.0 info: Table 'SEQUENCE' metadata: md5 ok
    # 2022-12-21T22:34:24 vdb-validate.2.10.0 info: Column 'ALTREAD': checksums ok
    # 2022-12-21T22:34:32 vdb-validate.2.10.0 info: Column 'QUALITY': checksums ok
    # 2022-12-21T22:34:37 vdb-validate.2.10.0 info: Column 'READ': checksums ok
    # 2022-12-21T22:34:37 vdb-validate.2.10.0 info: Column 'SPOT_GROUP': checksums ok
    # 2022-12-21T22:34:37 vdb-validate.2.10.0 info: Database 'SRR15048733/SRR15048733.sra' contains only unaligned reads
    # 2022-12-21T22:34:37 vdb-validate.2.10.0 info: Database 'SRR15048733.sra' is consistent

    # ... etc (similarly for all other accessions, unless of course there is a problem)

    # note: could send output to file if you didn't want to check it on screen

################################################################################################
# (6) Generate fastq files from prefetched sra files

    # N.B. if disk space is limiting, see Step (6*) below in which each sample's uncompressed
    # fastq files are gzipped before the next sample's fastq files are created

    ## this job takes a while and should be done in a batch job with a script something like this

    #    #!/bin/bash
    #    #SBATCH -N 1             # request one whole node (should have 24 cores)
    #    #SBATCH -n 6             # request 6 cores
    #    #SBATCH -t 0-12:00       # runtime in D-HH:MM
    #    #SBATCH -p shared        # partition to submit to
    #    #SBATCH -J fasterq-dump  # job name
    #
    #    source activate SRAtools
    #
    #    cd ~/permafrost-pathogens/data/fastq
    #
    #    for accession in $(<../../metadata/SraAccList.txt)
    #    do
    #        echo ${accession}
    #        fasterq-dump --threads 6 --temp . ${accession}/${accession}.sra
    #    done

################################################################################################
# (7) Compress fastq files in parallel on a compute node

    ## this job takes a while and should be done in a batch job with a script something like this

    ## not tested

    #    #!/bin/bash
    #    #SBATCH -N 1             # request one whole node (should have 24 cores)
    #    #SBATCH -t 0-12:00       # runtime in D-HH:MM
    #    #SBATCH -p shared        # partition to submit to
    #    #SBATCH -J gzip_fastqs   # job name
    #
    #    cd ~/permafrost-pathogens/data/fastq
    #
    #    # use gnu parallel to parallelize, since number of fastq files exceeds number of cores
    #    module load linuxbrew/colsa # to make gnu parallel available
    #    ls *.fastq | parallel 'gzip {}'

################################################################################################
# (8) Remove .sra files

    for accession in $(<SraAccList.txt)
    do
    rm ${accession}/${accession}.sra && rmdir ${accession}
    done

    # You should now be done with the data download.

################################################################################################
# (6*) Alternatively, in place of (6)-(8), loop through all the accessions to avoid having many
#      uncompressed fastq files at the end of step (6), since this may exhaust the disk quota.

    ## this job takes a while and should be done in a batch job with a script something like this

    # #!/bin/bash
    # #SBATCH -N 1             # request one whole node (should have 24 cores)
    # #SBATCH -n 6             # request six cores
    # #SBATCH --mem=16G        # request 16GB memory on the node
    # #SBATCH -t 0-12:00       # runtime in D-HH:MM
    # #SBATCH -p shared        # partition to submit to
    # #SBATCH -J fasterq-dump  # job name
    #
    # source activate SRAtools
    #
    # cd ~/permafrost-pathogens/data/fastq
    #
    # for accession in $(<../../metadata/SraAccList.txt)
    # do
    #     echo ${accession}
    #     # get fastq files; 6 threads is the default
    #         fasterq-dump --threads 6 --temp . ${accession}/${accession}.sra
    #     # compress fastq files
    #         gzip ${accession}.sra_1.fastq &
    #         gzip ${accession}.sra_2.fastq &
    #         wait
    #     # remove .sra files
    #         # rm ${accession}/${accession}.sra && rmdir ${accession}
    #    echo ${accession} complete
    # done

################################################################################################
# (9) optionally: calculate md5 hashes for downloaded (uncompressed) files to compare against saved values

    # Note: it may be better to rely on vdb-validate (see above)
    # -- it seems that the md5 hash for the uncompressed fastq data depends on the number of threads use
    # to run fasterq-dump, presumably because using a different number of threads has the effect of
    # reordering the fastq file (the default of 6 threads to generate the md5 sums at data/md5sums.txt)

    # Nonetheless, the code is included here in case it is useful.

    # Get an interactive compute session
        salloc -p shared -t 0-02:00
        # then connect to compute node using ssh

    conda activate SRAtools

    cd data/fastq

    # run in serial:

        rm -rf calculated_md5sums.txt # replace previous md5 sums if you have run this before
        for fastq in *.fastq.gz
        do
            echo -e "${fastq%.gz}\t"`zcat ${fastq} | md5sum | awk '{print $1}'` >>calculated_md5sums.txt
        done

    # or run in parallel:

        module load linuxbrew/colsa # to make gnu parallel available
        # calculate hashes
            ls *.fastq.gz | parallel "zcat {} | md5sum >{}.md5sum"
        # consolidate into a single file
            rm -rf calculated_md5sums.txt # replace previous md5 sums if you have run this before
            for md5 in *.md5sum
            do
                echo -e "${md5%.gz.md5sum}\t"`cat ${md5} | awk '{print $1}'` >>calculated_md5sums.txt
                rm ${md5}
            done

################################################################################################
