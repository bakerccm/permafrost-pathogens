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
import re
import os
METADATA = pd.read_csv(METADATA_FILE, sep = '\t', index_col = 'sample')
ALL_SAMPLES = list(METADATA.index)

# this is not a very clean way to do this:
#READ1_FILES = [re.sub(".fastq.gz", "", os.path.basename(file)) for file in METADATA['read1']]
#READ2_FILES = [re.sub(".fastq.gz", "", os.path.basename(file)) for file in METADATA['read2']]

################################
# default rules

# rule all:
#     input:
#         'out/multiqc'

#rule all_run_fastqc:
#    input:
#        ['out/fastqc/' + file + '_fastqc.zip' for file in READ1_FILES],
#        ['out/fastqc/' + file + '_fastqc.zip' for file in READ2_FILES]

################################

rule all_raw_data_links:
    input:
        expand('data/raw/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {'R1','R2'})

rule raw_data_link:
    input:
        read1 = lambda wildcards: RAW_DATA_DIR + "/" + METADATA.loc[wildcards.sample,'read1'],
        read2 = lambda wildcards: RAW_DATA_DIR + "/" + METADATA.loc[wildcards.sample,'read2']
    output:
        read1 = 'data/raw/{sample}_R1.fastq.gz',
        read2 = 'data/raw/{sample}_R2.fastq.gz'
    shell:
        '''
        ln -s {input.read1} {output.read1}
        ln -s {input.read2} {output.read2}
        '''

################################

# generate fastqc quality reports for each fastq.gz file

rule all_run_fastqc:
    expand('out/raw/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})

rule run_fastqc:
    input:
        'data/raw/{sample}_{read}.fastq.gz'
    output:
        'out/raw/{sample}_{read}_fastqc.html',
        'out/raw/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/raw {input}'

################################

# use multiQC to summarize fastqc results
# rule multiQC:
#     input:
#         fastqc = ['out/fastqc/' + file + '_fastqc.zip' for file in READ1_FILES] +
#                  ['out/fastqc/' + file + '_fastqc.zip' for file in READ2_FILES]
#         #star = expand('out/alignment/{sample}_Log.final.out', sample = list(SAMPLES.index)), # .Log.final.out files from STAR alignment
#         #featureCounts = 'out/featureCounts/bulkseq_featureCounts.txt.summary', # .summary file from featureCounts
#         #htseq_count = expand('out/htseq_count/{sample}.txt', sample = list(SAMPLES.index)) # .txt files from htseq-count
#     output:
#         directory('out/multiqc')
#     log:
#         'out/multiqc/multiqc_report.log'
#     conda:
#         'envs/multiqc.yaml'
#     shell:
#         'multiqc -o {output} {input} 2>{log}'

################################
