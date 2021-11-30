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
        'out/raw/multiqc'

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
# generate fastqc quality reports for each fastq.gz file

rule all_run_fastqc:
    input:
        expand('out/raw/fastqc/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})

rule run_fastqc:
    input:
        'data/links/{sample}_{read}.fastq.gz'
    output:
        'out/raw/fastqc/{sample}_{read}_fastqc.html',
        'out/raw/fastqc/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/raw/fastqc {input}'

################################
# use multiQC to summarize fastqc results
rule multiQC:
    input:
        fastqc = expand('out/raw/fastqc/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
        # star = expand('out/alignment/{sample}_Log.final.out', sample = list(SAMPLES.index)), # .Log.final.out files from STAR alignment
        # featureCounts = 'out/featureCounts/bulkseq_featureCounts.txt.summary', # .summary file from featureCounts
        # htseq_count = expand('out/htseq_count/{sample}.txt', sample = list(SAMPLES.index)) # .txt files from htseq-count
    output:
        directory('out/raw/multiqc')
    params:
        inputdir = 'out/raw/fastqc'
    log:
        'out/raw/multiqc/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc -o {output} {params.inputdir} 2>{log}'

################################
# remove sequencing adaptors

rule all_cutadapt:
    input:
        expand('out/cutadapt/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {'R1','R2'})

rule cutadapt:
    input:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    output:
        read1 = 'out/cutadapt/{sample}_R1.fastq.gz',
        read2 = 'out/cutadapt/{sample}_R2.fastq.gz'
    params:
        adapter_fwd = config['cutadapt']['adapter_fwd'],
        adapter_rev = config['cutadapt']['adapter_rev']
    conda:
        'envs/cutadapt-3.5.yaml'
    threads: 8
    shell:
        'cutadapt -a {params.adapter_fwd} -A {params.adapter_rev} -o {output.read1} -p {output.read2} {input.read1} {input.read2}'

# -n 2 -m {params.min_length}
################################


################################
