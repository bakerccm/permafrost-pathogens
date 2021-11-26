# snakefile for permafrost pathogens processing and analysis

# Chris Baker
# https://github.com/bakerccm
# 26 November 2021

################################
## get config file

configfile: 'config/config.yaml'
METADATA_FILE = config['sample_metadata']

################################
## wildcard constraints

wildcard_constraints:
    sample = '[^_/]+',
    fastq = '[^_/]+',
    file = '[^/\.]+'

################################
## get sample and reference genome metadata

import pandas as pd
METADATA = pd.read_csv(METADATA_FILE, sep = '\t', index_col = 'sample')
ALL_SAMPLES = list(METADATA.index)

################################
# default rules

rule all:
    input:
        'out/multiqc'

rule all_run_fastqc:
    input:
        ['out/fastqc/' + fastq + '_fastqc.zip' for fastq in METADATA['read1']],
        ['out/fastqc/' + fastq + '_fastqc.zip' for fastq in METADATA['read2']]

################################

# generate fastqc quality reports for each fastq.gz file
rule run_fastqc:
    input:
        'data/{fastq}.fastq.gz'
    output:
        'out/fastqc/{fastq}_fastqc.html',
        'out/fastqc/{fastq}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/fastqc {input}'

################################

# use multiQC to summarize fastqc results
rule multiQC:
    input:
        fastqc = ['out/fastqc/' + fastq + '_fastqc.zip' for fastq in METADATA['read1']] +
                 ['out/fastqc/' + fastq + '_fastqc.zip' for fastq in METADATA['read2']]
        #star = expand('out/alignment/{sample}_Log.final.out', sample = list(SAMPLES.index)), # .Log.final.out files from STAR alignment
        #featureCounts = 'out/featureCounts/bulkseq_featureCounts.txt.summary', # .summary file from featureCounts
        #htseq_count = expand('out/htseq_count/{sample}.txt', sample = list(SAMPLES.index)) # .txt files from htseq-count
    output:
        directory('out/multiqc')
    log:
        'out/multiqc/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc -o {output} {input} 2>{log}'

################################
