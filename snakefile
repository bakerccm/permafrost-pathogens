# snakefile for permafrost pathogens processing and analysis

# Chris Baker
# https://github.com/bakerccm
# 26 November 2021

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

################################
# default rules

rule all:
    input:
        'out/raw/multiqc_report.html', 'out/cutadapt/multiqc_report.html', 'out/bbduk/multiqc_report.html', 'out/sickle/multiqc_report.html', 'out/fastq-join/multiqc_report.html'

################################
rule all_raw_data_links:
    input:
        expand('data/raw/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {'R1','R2'})

rule raw_data_link:
    input:
        read1 = lambda wildcards: RAW_DATA_DIR + "/" + METADATA.loc[wildcards.sample,'read1'],
        read2 = lambda wildcards: RAW_DATA_DIR + "/" + METADATA.loc[wildcards.sample,'read2']
    output:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    shell:
        '''
        ln -s {input.read1} {output.read1}
        ln -s {input.read2} {output.read2}
        '''

################################
## QC for raw data files

# generate fastqc quality reports for each fastq.gz file

rule all_raw_fastqc:
    input:
        expand('out/raw/fastqc/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})

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

# use multiQC to summarize fastqc results
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
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
## use cutadapt to remove sequencing adaptors

rule all_cutadapt:
    input:
        expand('out/cutadapt/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {'R1','R2'})

# consider migrating to cutadapt wrapper:
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/cutadapt/pe.html
rule cutadapt:
    input:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    output:
        read1 = 'out/cutadapt/{sample}_R1.fastq.gz',
        read2 = 'out/cutadapt/{sample}_R2.fastq.gz',
        qc = 'out/cutadapt/{sample}.qc.txt'
    params:
        quality_cutoff = config['cutadapt']['quality_cutoff'],
        adapter_fwd = config['cutadapt']['adapter_fwd'],
        adapter_rev = config['cutadapt']['adapter_rev'],
        max_error_rate = config['cutadapt']['max_error_rate'],
        minimum_length = config['cutadapt']['minimum_length']
    conda:
        'envs/cutadapt-3.5.yaml'
    threads: 4
    shell:
        '''
        cutadapt -a {params.adapter_fwd} -A {params.adapter_rev} \
        -j {threads} -e {params.max_error_rate} -m {params.minimum_length} -q {params.quality_cutoff}\
        -o {output.read1} -p {output.read2} {input.read1} {input.read2} >{output.qc}
        '''

rule cutadapt_fastqc:
    input:
        'out/cutadapt/{sample}_{read}.fastq.gz'
    output:
        'out/cutadapt/{sample}_{read}_fastqc.html',
        'out/cutadapt/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/cutadapt {input}'

rule cutadapt_multiQC:
    input:
        fastqc = expand('out/cutadapt/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'}),
        cutadapt = expand('out/cutadapt/{sample}.qc.txt', sample = ALL_SAMPLES)
    output:
        'out/cutadapt/multiqc_report.html'
    params:
        inputdir = 'out/cutadapt',
        outputdir = 'out/cutadapt'
    log:
        'out/cutadapt/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
rule bbduk:
    input:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    output:
        read1 = 'out/bbduk/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk/{sample}_R2.fastq.gz'
    params:
        memory = config['bbduk']['memory'],
        ref = config['bbduk']['ref'],
        ktrim = config['bbduk']['ktrim'],
        k = config['bbduk']['k'],
        mink = config['bbduk']['mink'],
        hdist = config['bbduk']['hdist'],
        trim_params = config['bbduk']['trim_params'],
        qtrim = config['bbduk']['qtrim'],
        trimq = config['bbduk']['trimq'],
        minlength = config['bbduk']['minlength']
    log:
        'out/bbduk/{sample}.log'
    conda:
        'envs/bbtools.yaml'
    shell:
        '''
        bbduk.sh {params.memory} in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} \
        ref={params.ref} ktrim={params.ktrim} k={params.k} {params.trim_params} \
        qtrim={params.qtrim} trimq={params.trimq} minlength={params.minlength} \
        &>>{log}
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
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
## use sickle to trim and filter reads

# consider migrating to sickle wrapper:
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/sickle/pe.html
rule sickle:
    input:
        read1 = 'out/cutadapt/{sample}_R1.fastq.gz',
        read2 = 'out/cutadapt/{sample}_R2.fastq.gz'
    output:
        read1 = 'out/sickle/{sample}_R1.fastq.gz',
        read2 = 'out/sickle/{sample}_R2.fastq.gz',
        unpaired = 'out/sickle/{sample}_unpaired.fastq.gz'
    log:
        'out/sickle/{sample}_sickle.log'
    params:
        quality_threshold = config['sickle']['quality_threshold'],
        length_threshold = config['sickle']['length_threshold']
    conda:
        'envs/sickle.yaml'
    threads: 4
    shell:
        '''
        sickle pe \
            -f {input.read1} -r {input.read2} \
            -t sanger -g \
            -o {output.read1} -p {output.read2} -s {output.unpaired} \
            -q {params.quality_threshold} -l {params.length_threshold} \
            > {log}
        '''

rule sickle_fastqc:
    input:
        'out/sickle/{sample}_{read}.fastq.gz'
    output:
        'out/sickle/{sample}_{read}_fastqc.html',
        'out/sickle/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/sickle {input}'

rule sickle_multiQC:
    input:
        fastqc = expand('out/sickle/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'}),
        sickle = expand("out/sickle/{sample}_sickle.log", sample = ALL_SAMPLES)
    output:
        'out/sickle/multiqc_report.html'
    params:
        inputdir = 'out/sickle',
        outputdir = 'out/sickle'
    log:
        'out/sickle/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
## use fastq-join from ea-utils to merge paired end reads
# see https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqJoin.md

rule fastq_join:
    input:
        read1 = 'out/bbduk/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk/{sample}_R2.fastq.gz'
    output:
        join = 'out/fastq-join/{sample}_join.fastq.gz',
        unmatched1 = 'out/fastq-join/{sample}_un1.fastq.gz',
        unmatched2 = 'out/fastq-join/{sample}_un2.fastq.gz'
    params:
        max_percent_difference = config['fastq-join']['max_percent_difference'],
        min_overlap = config['fastq-join']['min_overlap']
    conda:
        'envs/ea-utils.yaml'
    shell:
        '''
        fastq-join -p {params.max_percent_difference} -m {params.min_overlap} \
        {input.read1} {input.read2} \
        -o {output.unmatched1} -o {output.unmatched2} -o {output.join}
        '''

rule fastq_join_fastqc:
    input:
        'out/fastq-join/{sample}_join.fastq.gz'
    output:
        'out/fastq-join/{sample}_join_fastqc.html',
        'out/fastq-join/{sample}_join_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/fastq-join {input}'

rule fastq_join_multiQC:
    input:
        fastqc = expand('out/fastq-join/{sample}_join_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
    output:
        'out/fastq-join/multiqc_report.html'
    params:
        inputdir = 'out/fastq-join',
        outputdir = 'out/fastq-join'
    log:
        'out/fastq-join/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# megahit assembly

rule megahit:
    input:
        read1 = ["out/bbduk/35m-t0-R1_R1.fastq.gz","out/bbduk/35m-t0-R2_R1.fastq.gz"],
        read2 = ["out/bbduk/35m-t0-R1_R2.fastq.gz","out/bbduk/35m-t0-R2_R2.fastq.gz"]
    output:
        directory("out/megahit")
    threads: 16
    conda:
        'envs/megahit.yaml'
    shell:
        '''
        # get file names and save as array
            read1_files_space=({input.read1})
            read2_files_space=({input.read2})
        
        # use sed to replace spaces with commas (assumes no spaces in file names)
            read1_files_comma=`echo ${read1_files_space[@]} | sed 's/ /,/g'`
            read2_files_comma=`echo ${read2_files_space[@]} | sed 's/ /,/g'`
        
        # uses printf to join array with commas (works even if file names contain spaces, but these will be passed through and will mess up snakemake unless file names are quoted)
            # printf -v joined_read1 '%s,' "${read1_files_space[@]}"
            # printf -v joined_read2 '%s,' "${read2_files_space[@]}"
            # read1_files_comma=`echo "${joined_read1%,}"`
            # read2_files_comma=`echo "${joined_read2%,}"`
        
        megahit -1 ${read1_files_comma} -2 ${read2_files_comma} -t {threads} -o {output}
        '''

################################
# metaspades assembly

rule metaspades:
    input:
        # may need to use yamle file or some other approach if co-assmblying many files
        pe1_1 = "out/bbduk/35m-t0-R1_R1.fastq.gz", # library 1, read 1
        pe1_2 = "out/bbduk/35m-t0-R1_R2.fastq.gz", # library 1, read 2
        pe2_1 = "out/bbduk/35m-t0-R2_R1.fastq.gz", # library 1, read 1
        pe2_2 = "out/bbduk/35m-t0-R2_R2.fastq.gz"  # library 2, read 2
    output:
        directory("out/metaspades")
    threads: 16
    conda:
        'envs/spades.yaml'
    shell:
        '''
        metaspades.py -t {threads} \
        --pe1-1 {input.pe1_1} --pe1-2 {input.pe1_2} \
        --pe2-1 {input.pe2_1} --pe2-2 {input.pe2_2} \
        -o {output}
        '''

################################
