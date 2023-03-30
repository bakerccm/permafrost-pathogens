# snakefile for permafrost pathogens processing and analysis

# Chris Baker
# https://github.com/bakerccm
# 29 December 2022

from snakemake.utils import min_version
min_version("6.4.1")

################################
## get config file

configfile: 'config/config.yaml'
METADATA_FILE = config['sample_metadata']
RAW_DATA_DIR = config['raw_data_dir']

################################
## wildcard constraints

wildcard_constraints:
    sample = '[^_/]+',
    file = '[^/]+'

################################
## get sample and reference genome metadata

import pandas as pd
METADATA = pd.read_csv(METADATA_FILE, sep = '\t', index_col = 'sample')
ALL_SAMPLES = list(METADATA.index)
GOOD_SAMPLES = list(METADATA[METADATA.co_assembly != "NA"].index) # samples shown as "." in the co_assembly column are controls or poor quality --> exclude from assemblies

################################
# default rule

rule all:
    input:
        'out/raw/multiqc_report.html',
        'out/singlem/all_good_samples.otu_table.tsv',
        'out/bbduk/multiqc_report.html',
        'out/bbduk_noPhiX/multiqc_report.html',
        'out/bbduk_noPhiX_fastuniq/multiqc_report.html',
        'out/phyloflash/all_good_samples.phyloFlash_compare.barplot.pdf',
        expand('out/metaquast/{assembly}/report.txt', assembly = {'35m','45m', '60m', '83m', 'NT'}),
        expand('out/maxbin2/{assembly}/done', assembly = {'35m','45m', '60m', '83m', 'NT'}),
        expand('out/maxbin2_checkm/{assembly}.txt', assembly = {'35m','45m', '60m', '83m', 'NT'}),
        # expand these to all assemblies and bins later
            'out/maxbin2_prokka/35m/35m.001',
            'out/maxbin2_staramr/35m/35m.001',
            'out/maxbin2_rgi/35m/35m.001.txt'

################################
## make links to raw data files (note: renames files to reflect sample labels)

rule all_raw_data_links:
    input:
        expand('data/links/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {'R1','R2'})

rule raw_data_link:
    input:
        # note these are specified relative to snakemake root
        read1 = lambda wildcards: RAW_DATA_DIR + "/" + METADATA.loc[wildcards.sample,'read1'],
        read2 = lambda wildcards: RAW_DATA_DIR + "/" + METADATA.loc[wildcards.sample,'read2']
    output:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    shell:
        # note use of ../../ in target to ensure correct location relative to link
        '''
        ln -s ../../{input.read1} {output.read1}
        ln -s ../../{input.read2} {output.read2}
        '''

################################
## qc for raw data files

rule all_raw_qc:
    input:
        'out/raw/multiqc_report.html'

rule raw_fastqc:
    input:
        'data/links/{sample}_{read}.fastq.gz'
    output:
        'out/raw/{sample}_{read}_fastqc.html',
        'out/raw/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/raw {input}'

rule raw_multiQC:
    input:
        fastqc = expand('out/raw/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
    output:
        'out/raw/multiqc_report.html'
    params:
        inputdir = 'out/raw',
        outputdir = 'out/raw'
    log:
        'out/raw/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -f -o {params.outputdir} {params.inputdir} 2>{log}'

################################
## run rules bbduk, bbduk_noPhiX and fastuniq, plus QC for each step

# note: these rules should be run together like this if the outputs of bbduk and bbduk_noPhiX are temp()

rule bbduk_fastuniq_all:
    input:
        "out/bbduk/multiqc_report.html",
        "out/bbduk_noPhiX/multiqc_report.html",
        "out/bbduk_noPhiX_fastuniq/multiqc_report.html"

################################
## trim adapters etc using bbduk

# see bbduk overview here http://seqanswers.com/forums/showthread.php?t=42776
rule bbduk:
    input:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    output:
        # consider making these temp() but only if it is possible to run bbduk,
        # bbduk_noPhiX and fastuniq (plus QC) together in a single job
        #    read1 = temp('out/bbduk/{sample}_R1.fastq.gz'),
        #    read2 = temp('out/bbduk/{sample}_R2.fastq.gz'),
        read1 = 'out/bbduk/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk/{sample}_R2.fastq.gz',
        stats = 'out/bbduk/{sample}_stats.txt'
    params:
        memory = config['bbduk']['memory'],
        ref = config['bbduk']['adapter_trimming']['ref'],
        ktrim = config['bbduk']['adapter_trimming']['ktrim'],
        k = config['bbduk']['adapter_trimming']['k'],
        mink = config['bbduk']['adapter_trimming']['mink'],
        hdist = config['bbduk']['adapter_trimming']['hdist'],
        trim_params = config['bbduk']['adapter_trimming']['trim_params'],
        qtrim = config['bbduk']['quality_filtering']['qtrim'],
        trimq = config['bbduk']['quality_filtering']['trimq'],
        minlength = config['bbduk']['quality_filtering']['minlength']
    log:
        'out/bbduk/{sample}.log'
    conda:
        'envs/bbtools.yaml'
    threads:
        config['bbduk']['threads']
    shell:
        '''
        bbduk.sh {params.memory} threads={threads} \
        in1={input.read1} in2={input.read2} \
        out1={output.read1} out2={output.read2} \
        ref={params.ref} ktrim={params.ktrim} k={params.k} mink={params.mink} hdist={params.hdist} {params.trim_params} \
        qtrim={params.qtrim} trimq={params.trimq} minlength={params.minlength} \
        stats={output.stats} &>>{log}
        '''

rule bbduk_fastqc:
    input:
        'out/bbduk/{sample}_{read}.fastq.gz'
    output:
        'out/bbduk/{sample}_{read}_fastqc.html',
        'out/bbduk/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/bbduk {input}'

rule bbduk_multiQC:
    input:
        fastqc = expand('out/bbduk/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
    output:
        'out/bbduk/multiqc_report.html'
    params:
        inputdir = 'out/bbduk',
        outputdir = 'out/bbduk'
    log:
        'out/bbduk/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -f -o {params.outputdir} {params.inputdir} 2>{log}'

################################
## remove any residual PhiX using bbduk

rule bbduk_noPhiX:
    input:
        read1 = 'out/bbduk/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk/{sample}_R2.fastq.gz'
    output:
        # consider making these temp() but only if it is possible to run bbduk,
        # bbduk_noPhiX and fastuniq (plus QC) together in a single job
        #    read1_unmatched = temp('out/bbduk_noPhiX/{sample}_unmatched_R1.fastq.gz'), # unmatched reads are not PhiX
        #    read2_unmatched = temp('out/bbduk_noPhiX/{sample}_unmatched_R2.fastq.gz'),
        read1_unmatched = 'out/bbduk_noPhiX/{sample}_unmatched_R1.fastq.gz', # unmatched reads are not PhiX
        read2_unmatched = 'out/bbduk_noPhiX/{sample}_unmatched_R2.fastq.gz',
        stats = 'out/bbduk_noPhiX/{sample}_stats.txt'
    params:
        memory = config['bbduk']['memory'],
        ref = config['bbduk']['contaminant_filtering']['ref'],
        k = config['bbduk']['contaminant_filtering']['k'],
        hdist = config['bbduk']['contaminant_filtering']['hdist']
    log:
        'out/bbduk_noPhiX/{sample}.log'
    conda:
        'envs/bbtools.yaml'
    threads:
        config['bbduk']['threads']
    shell:
        '''
        bbduk.sh {params.memory} threads={threads} \
        in1={input.read1} in2={input.read2} \
        out1={output.read1_unmatched} out2={output.read2_unmatched} \
        ref={params.ref} k={params.k} hdist={params.hdist} \
        stats={output.stats} &>>{log}
        '''

rule bbduk_noPhiX_fastqc:
    input:
        'out/bbduk_noPhiX/{sample}_unmatched_{read}.fastq.gz'
    output:
        'out/bbduk_noPhiX/{sample}_unmatched_{read}_fastqc.html',
        'out/bbduk_noPhiX/{sample}_unmatched_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/bbduk_noPhiX {input}'

rule bbduk_noPhiX_multiQC:
    input:
        fastqc = expand('out/bbduk_noPhiX/{sample}_unmatched_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
    output:
        'out/bbduk_noPhiX/multiqc_report.html'
    params:
        inputdir = 'out/bbduk_noPhiX',
        outputdir = 'out/bbduk_noPhiX'
    log:
        'out/bbduk_noPhiX/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -f -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# deduplicate using FastUniq

# commands for testing
# snakemake --use-conda --rerun-incomplete -j 16 out/bbduk_noPhiX_fastuniq/35m-t0-R2_{R1,R2}.fastq.gz
# snakemake --use-conda --rerun-incomplete -j 56 -np fastuniq_all

rule fastuniq_decompress_inputs:
    input:
        'out/bbduk_noPhiX/{sample}_unmatched_{read}.fastq.gz' # unmatched reads are not PhiX
    output:
        temp('out/bbduk_noPhiX/{sample}_unmatched_{read}.fastq')
    shell:
        'zcat {input} >{output}'

rule fastuniq:
    input:
        read1 = 'out/bbduk_noPhiX/{sample}_unmatched_R1.fastq', # unmatched reads are not PhiX
        read2 = 'out/bbduk_noPhiX/{sample}_unmatched_R2.fastq'  # note inputs are uncompressed
    output:
        filelist = 'out/bbduk_noPhiX_fastuniq/{sample}_input_filelist.txt',
        read1 = 'out/bbduk_noPhiX_fastuniq/{sample}_R1.fastq', # note outputs are uncompressed
        read2 = 'out/bbduk_noPhiX_fastuniq/{sample}_R2.fastq'  # but these are actually removed if rule fastuniq_compress_outputs is run
    params:
        filelist_name = 'out/bbduk_noPhiX_fastuniq/{sample}_input_filelist.txt'
    conda:
        'envs/fastuniq.yaml'
    threads:
        config['fastuniq']['threads']
    shell:
        '''
        # save input file names to file
            echo {input.read1} >{output.filelist}
            echo {input.read2} >>{output.filelist}
        # run fastuniq
            fastuniq -i {output.filelist} -t q -c 0 \
            -o {output.read1} -p {output.read2}
        '''

rule fastuniq_compress_outputs:
    input:
        'out/bbduk_noPhiX_fastuniq/{sample}_{read}.fastq'
    output:
        'out/bbduk_noPhiX_fastuniq/{sample}_{read}.fastq.gz'
    shell:
        'gzip {input}'

rule fastuniq_fastqc:
    input:
        'out/bbduk_noPhiX_fastuniq/{sample}_{read}.fastq.gz'
    output:
        'out/bbduk_noPhiX_fastuniq/{sample}_{read}_fastqc.html',
        'out/bbduk_noPhiX_fastuniq/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/bbduk_noPhiX_fastuniq {input}'

rule fastuniq_multiQC:
    input:
        fastqc = expand('out/bbduk_noPhiX_fastuniq/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
    output:
        'out/bbduk_noPhiX_fastuniq/multiqc_report.html'
    params:
        inputdir = 'out/bbduk_noPhiX_fastuniq',
        outputdir = 'out/bbduk_noPhiX_fastuniq'
    log:
        'out/bbduk_noPhiX_fastuniq/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -f -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# megahit co-assembly

# co-assembles samples excluding failed samples, according to co_assembly column in samples.tsv metadata file
# runtimes:
#   35m: 04:02:15
#   45m: 07:22:25
#   60m: 08:35:44
#   83m: 05:29:18
#   NT : 10:20:45
rule megahit_coassembly:
    input:
        read1 = lambda wildcards: ["out/bbduk_noPhiX_fastuniq/" + sample + "_R1.fastq.gz" for sample in list(METADATA[METADATA.co_assembly == wildcards.assembly].index)],
        read2 = lambda wildcards: ["out/bbduk_noPhiX_fastuniq/" + sample + "_R2.fastq.gz" for sample in list(METADATA[METADATA.co_assembly == wildcards.assembly].index)]
    output:
        "out/megahit/{assembly}/final.contigs.fa"
    params:
        output_dir = "out/megahit/{assembly}"
    log:
        "out/megahit/{assembly}.log"
    threads:
        config['megahit']['threads']
    conda:
        'envs/megahit.yaml'
    shell:
        '''
        # get file names and store as array
            read1_files_space=({input.read1})
            read2_files_space=({input.read2})

        # use sed to replace spaces with commas (assumes no spaces in file names)
            read1_files_comma=`echo ${{read1_files_space[@]}} | sed 's/ /,/g'`
            read2_files_comma=`echo ${{read2_files_space[@]}} | sed 's/ /,/g'`

        # write files names to log file
            echo "Co-assembly: {wildcards.assembly}" >{log}
            echo "Read 1 files" >>{log}
            echo $read1_files_comma >>{log}
            echo "Read 2 files" >>{log}
            echo $read2_files_comma >>{log}

        # remove output directory (megahit will fail if already exists)
            rm -rf {params.output_dir}

        megahit -1 ${{read1_files_comma}} -2 ${{read2_files_comma}} -t {threads} -o {params.output_dir}
        '''

# takes 5 min to run without gene finding
rule megahit_metaquast:
    input:
        "out/megahit/{assembly}/final.contigs.fa"
    output:
        "out/metaquast/{assembly}/report.txt" # there are several other output files in this directory too
    params:
        output_dir = "out/metaquast/{assembly}"
    threads:
        config['metaquast']['threads']
    conda:
        'envs/quast.yaml'
    shell:
        # use --gene-finding (of -f) to find genes using MetaGeneMark
        # -- but this requires either a sourceforge installation of quast, or a manual change to the bioconda installation
        # -- see https://github.com/ablab/quast/issues/84
        # use --max-ref-number 0 to skip searching against SILVA and downloading refs
        'metaquast.py -o {params.output_dir} -t {threads} --max-ref-number 0 {input}'

################################
# this rule fails with errors messages about missing contextvars and missing django

# rule reformat_megahit_contigs:
#     input:
#         "out/megahit/{assembly}/final.contigs.fa"
#     output:
#         "out/megahit/{assembly}/contigs.fa"
#     params:
#         min_length = 1000
#     conda:
#         'envs/anvio-minimal.yaml'
#     shell:
#         'anvi-script-reformat-fasta {input} -o {output} -l {params.min_length} --simplify-names'

################################
# use phyloflash to do some SSU-based analysis

# Note: phyloflash database should be generated first following instructions at
#     http://hrgv.github.io/phyloFlash/install.html

# Ideally, you can use phyloFlash_makedb.pl to automate the process:
#    # activate conda environment with phyloflash in it
#    # then ...
#    mkdir -p databases/phyloflash
#    cd databases/phyloflash
#    phyloFlash_makedb.pl --remote

# This should automatically retrieve sequence files from the internet and build databases.

# Sometimes it fails though, in which case the files can be downloaded manually and supplied
# to phyloFlash_makedb.pl (follow instructions at http://hrgv.github.io/phyloFlash/install.html).

# Note that the location of the database (e.g. databases/phyloflash/138.1) needs to be 
# supplied in /config/config.yaml

rule all_phyloflash:
    input:
        expand('out/phyloflash/{sample}.phyloFlash.html', sample = GOOD_SAMPLES)

rule phyloflash:
    input:
        read1 = 'out/bbduk_noPhiX_fastuniq/{sample}_R1.fastq.gz', # cleaned data files
        read2 = 'out/bbduk_noPhiX_fastuniq/{sample}_R2.fastq.gz'
    output:
        'out/phyloflash/{sample}.phyloFlash.html',
        'out/phyloflash/{sample}.phyloFlash.log',
        'out/phyloflash/{sample}.phyloFlash.tar.gz'
    params:
        output_dir = 'out/phyloflash',
        readlength = config['phyloflash']['readlength'],
        dbhome_dir = config['phyloflash']['dbhome_dir']
    threads:
        config['phyloflash']['threads']
    conda:
        'envs/phyloflash.yaml'
    shell:
        '''
        cd {params.output_dir}
        phyloFlash.pl -lib {wildcards.sample} -read1 ../../{input.read1} -read2 ../../{input.read2} \
        -CPUs {threads} -readlength {params.readlength} -dbhome ../../{params.dbhome_dir} \
        -emirge -poscov -treemap -zip -log
        '''

rule phyloflash_compare:
    input:
        expand('out/phyloflash/{sample}.phyloFlash.tar.gz', sample = GOOD_SAMPLES)
    output:
        'out/phyloflash/all_good_samples.phyloFlash_compare.barplot.pdf',
        'out/phyloflash/all_good_samples.phyloFlash_compare.heatmap.pdf',
        'out/phyloflash/all_good_samples.phyloFlash_compare.matrix.tsv',
        'out/phyloflash/all_good_samples.phyloFlash_compare.ntu_table.tsv'
    conda:
        'envs/phyloflash.yaml'
    shell:
        '''
        cd out/phyloflash
        phyloFlash_compare.pl --allzip --task barplot,heatmap,matrix,ntu_table -out all_good_samples.phyloFlash_compare
        # (or specify the actual files instead of the directory)
        '''

################################

rule all_singlem:
    input:
        expand('out/singlem/{sample}.otu_table.tsv', sample = GOOD_SAMPLES)

rule singlem:
    input:
        read1 = 'data/links/{sample}_R1.fastq.gz', # may want to use adapter-filtered reads instead; manual suggests not quality filtering since it can make the reads too short
        read2 = 'data/links/{sample}_R2.fastq.gz'
    output:
        'out/singlem/{sample}.otu_table.tsv'
    conda:
        'envs/singlem.yaml'
    threads:
        config['singlem']['pipe']['threads']
    shell:
        '''
        singlem pipe --forward {input.read1} --reverse {input.read2} \
        --otu_table {output} --threads {threads}
        '''
# other popular options:
# --output_extras	Output more detailed information in the OTU table.
# --assignment_method {pplacer,diamond,diamond_example}
#   Specify taxonomic assignment method [default: pplacer].
# conda activate /work2/08186/cbaker/frontera/permafrost-pathogens/.snakemake/conda/89cb197209c9f30f798aec2bc588a442
# see singlem pipe --full_help for more help

rule singlem_summarise:
    input:
        expand('out/singlem/{sample}.otu_table.tsv', sample = GOOD_SAMPLES)
    output:
        krona = 'out/singlem/all_good_samples.krona.html',
        OTU_table = 'out/singlem/all_good_samples.otu_table.tsv'
    conda:
        'envs/singlem.yaml'
    shell:
        '''
        singlem summarise --input_otu_tables {input} --krona {output.krona}
        singlem summarise --input_otu_tables {input} --output_otu_table {output.OTU_table}
        '''

# consider clustering OTUs, rarefying, beta diversity etc
################################
# create bowtie2 index for each megahit co-assembly

rule all_bowtie2_build:
    input:
        expand('out/megahit/{assembly}/bowtie2_index.1.bt2', assembly = {'35m','45m','60m','83m','NT'}) # just one of the files in each index

rule bowtie2_build:
    input:
        'out/megahit/{assembly}/final.contigs.fa'
    output:
        'out/megahit/{assembly}/bowtie2_index.1.bt2' # just one of the files; is this a good naming scheme?
    params:
        bt2_index = 'out/megahit/{assembly}/bowtie2_index' # file path stem for bowtie2 index
    conda:
        'envs/bowtie2.yaml'
    shell:
        'bowtie2-build {input} {params.bt2_index}'

rule all_bowtie2_mapping:
    input:
        expand('out/megahit/{assembly}/bowtie2_mapping/{sample}.bam', zip, assembly = list(METADATA.co_assembly[GOOD_SAMPLES]), sample = GOOD_SAMPLES)

rule bowtie2_mapping:
    input:
        read1 = 'out/bbduk_noPhiX_fastuniq/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk_noPhiX_fastuniq/{sample}_R2.fastq.gz',
        bowtie2_index = 'out/megahit/{assembly}/bowtie2_index.1.bt2' # just one of the files
    output:
        temp('out/megahit/{assembly}/bowtie2_mapping/{sample}.sam')
    log:
        'out/megahit/{assembly}/bowtie2_mapping/{sample}.log'
    params:
        bt2_index = 'out/megahit/{assembly}/bowtie2_index' # file path stem for bowtie2 index
    conda:
        'envs/bowtie2.yaml'
    threads:
        config['bowtie2']['threads']
    shell:
        '''
        bowtie2 -1 {input.read1} -2 {input.read2} -q \
        -x {params.bt2_index} --no-unal --threads {threads} -S {output} 2>{log}
        '''

rule bowtie2_compress_sort_index:
    input:
        'out/megahit/{assembly}/bowtie2_mapping/{sample}.sam'
    output:
        'out/megahit/{assembly}/bowtie2_mapping/{sample}.bam'
    conda:
        'envs/samtools.yaml'
    params:
        temp_unsorted_bam = 'out/megahit/{assembly}/bowtie2_mapping/{sample}_raw.bam'
    shell:
        '''
        # convert to unsorted bam
            samtools view -b -o {params.temp_unsorted_bam} {input}
        # sort bam file
            samtools sort -o {output} {params.temp_unsorted_bam}
            rm {params.temp_unsorted_bam}
        # index sorted bam
            samtools index {output}
        '''

################################

# get abundance info for each bowtie mapping

rule all_bowtie_coverage:
    input:
        expand('out/megahit/{assembly}/bowtie2_mapping/{sample}_coverage.txt', zip, assembly = list(METADATA.co_assembly[GOOD_SAMPLES]), sample = GOOD_SAMPLES),
        expand('out/megahit/{assembly}/bowtie2_mapping/{sample}_meandepth.txt', zip, assembly = list(METADATA.co_assembly[GOOD_SAMPLES]), sample = GOOD_SAMPLES)

rule get_bowtie_coverage:
    input:
        'out/megahit/{assembly}/bowtie2_mapping/{sample}.bam'
    output:
        coverage = 'out/megahit/{assembly}/bowtie2_mapping/{sample}_coverage.txt', # full samtools coverage output
        meandepth = 'out/megahit/{assembly}/bowtie2_mapping/{sample}_meandepth.txt' # just the meandepth column from the coverage output table
    conda:
        'envs/samtools.yaml'
    params:
        bamfile_list = 'out/megahit/{assembly}/bowtie2_mapping/{sample}_bamfile_list.txt'
    shell:
        '''
        echo {input} >{params.bamfile_list}
        samtools coverage -b {params.bamfile_list} > {output.coverage}
        awk '{{print $1"\t"$6}}' {output.coverage} | grep -v '^#' > {output.meandepth}
        '''

rule maxbin2:
    input:
        meandepth_tables = lambda wildcards: ('out/megahit/' + wildcards.assembly + '/bowtie2_mapping/' + sample + '_meandepth.txt' for sample in list(METADATA[METADATA["co_assembly"] == wildcards.assembly].index)),
        contigs = 'out/megahit/{assembly}/final.contigs.fa'
    output:
        'out/maxbin2/{assembly}/done'
    params:
        meandepth_filelist = 'out/maxbin2/{assembly}/{assembly}_meandepth_filelist.txt',
        output_file_header = 'out/maxbin2/{assembly}/{assembly}'
    conda:
        'envs/maxbin2.yaml'
    threads:
        config['maxbin2']['threads']
    shell:
        '''
        # get file names and store as array
            meandepth_tables=({input.meandepth_tables})

        # print file names to file separated by newlines
            printf "%s\n" "${{meandepth_tables[@]}}" > {params.meandepth_filelist}

        run_MaxBin.pl -contig {input.contigs} -out {params.output_file_header} \
        -thread {threads} -abund_list {params.meandepth_filelist}
        
        touch {output}
        '''

# --- at least one of the following parameters is needed
# (semi-required) -abund (contig abundance files. To be explained in Abundance session below.)
# (semi-required) -reads (reads file in fasta or fastq format. To be explained in Abundance session below.)
# (semi-required) -abund_list (a list file of all contig abundance files.)
# (semi-required) -reads_list (a list file of all reads file.)
# other options
# (optional) -prob_threshold (minimum probability for EM algorithm; default 0.8)
# (optional) -plotmarker (specify this option if you want to plot the markers in each contig. Installing R is a must for this option to work.)
# (optional) -verbose (as is. Warning: output log will be LOOOONG.)
# (optional) -markerset (choose between 107 marker genes by default or 40 marker genes. see Marker Gene Note for more information.)

# use checkM to assess putative genomes
# note databases need to be downloaded first. see install instructions https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm
#   mkdir -p databases/checkm
#   cd databases/checkm
#   wget 'https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz' # 276 MB download
#   tar -xvf 'checkm_data_2015_01_16.tar.gz' # unpacks to a 1.4G directory
#   rm 'checkm_data_2015_01_16.tar.gz'

# it would be good to make checkm output this to a file ... -f flag should do it I think

rule checkm:
    input:
        'out/maxbin2/{assembly}/done'
    output:
        dir = directory('out/maxbin2_checkm/{assembly}'), # checkm output folder
        file = 'out/maxbin2_checkm/{assembly}.txt'
    params:
        input_dir = 'out/maxbin2/{assembly}',
        database_dir = 'databases/checkm'
    conda:
        'envs/checkm.yaml'
    threads:
        config['checkm']['threads']
    shell:
        '''
        # set database directory
            checkm data setRoot {params.database_dir}
        # run checkm
        # note fasta extension is specified here - updated if using a different binning program
            checkm lineage_wf -t {threads} -x fasta {params.input_dir} -f {output.file} {output.dir}
        '''

# checkm lineage_wf runs the four mandatory steps of the lineage-specific workflow:
#   (M) > checkm tree <bin folder> <output folder>
#   (R) > checkm tree_qa <output folder>
#   (M) > checkm lineage_set <output folder> <marker file>
#   (M) > checkm analyze <marker file> <bin folder> <output folder>
#   (M) > checkm qa <marker file> <output folder>

# possible additional analyses see https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree
#   checkm tree_qa <tree folder>
# and possibly some additional qa outputs?
# see https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree
################################

# need to automatically work out what files to ask for here, but for now this will do:
# snakemake -j 56 --use-conda out/maxbin2_prokka/35m/35m.{001..053} # done
# snakemake -j 56 --use-conda out/maxbin2_prokka/45m/45m.{001..107} # done
# snakemake -j 56 --use-conda out/maxbin2_prokka/60m/60m.{001..099} # done
# snakemake -j 56 --use-conda out/maxbin2_prokka/83m/83m.{001..071} # done
# snakemake -j 56 --use-conda out/maxbin2_prokka/NT/NT.{001..096} # working on it
rule prokka:
    input:
        # note sure how to specify fasta(s) as output of maxbin2 rule since number of bins is not predefined
        # 'out/maxbin2/{assembly}/{bin}.fasta'
        'out/maxbin2/{assembly}/done'
    output:
        directory('out/maxbin2_prokka/{assembly}/{bin}') # maybe alter this once we know what output files look like
    params:
        input_fasta = 'out/maxbin2/{assembly}/{bin}.fasta'
        # ?? out = 'databases/checkm'
    conda:
        'envs/prokka.yaml'
    shell:
        '''
        prokka --outdir {output} --prefix {wildcards.bin} {params.input_fasta}
        '''

################################

rule staramr:
    input:
        'out/maxbin2/{assembly}/done'
    output:
        directory('out/maxbin2_staramr/{assembly}/{bin}')
    params:
        input_fasta = 'out/maxbin2/{assembly}/{bin}.fasta'
    conda:
        'envs/staramr.yaml'
    shell:
        '''
        staramr search -o {output} {params.input_fasta}
        '''

################################
# database needs to be set up first -- see https://github.com/arpcard/rgi#rgi-usage-documentation
#   mkdir -p databases/RGI
#   cd databases/RGI
#   wget https://card.mcmaster.ca/latest/data
#   tar -xvf data ./card.json
#   cd ../..

# rgi load --card_json /path/to/card.json --local

rule rgi:
    input:
        'out/maxbin2/{assembly}/done'
    output:
        'out/maxbin2_rgi/{assembly}/{bin}.txt'
    params:
        input_fasta = 'out/maxbin2/{assembly}/{bin}.fasta'
    conda:
        'envs/rgi.yaml'
    threads:
        config['rgi']['threads']
    shell:
        '''
        rgi main --input_sequence {params.input_fasta} \
          --output_file {output} --input_type contig --local \
          --low_quality --include_loose --clean --num_threads {threads} --split_prodigal_jobs
        '''

################################
