# snakefile template

# Chris Baker
# https://github.com/bakerccm
# DD MMMM YYYY

################################
## get config file

configfile: 'config/config.yaml'
METADATA_FILE = config['sample_metadata']

################################
## wildcard constraints

wildcard_constraints:
    sample = '[^_/]+',
    file = '[^/\.]+'

################################
## get sample and reference genome metadata

import pandas as pd
METADATA = pd.read_csv(METADATA_FILE, sep = '\t', index_col = 'sample')
ALL_SAMPLES = list(METADATA[(METADATA['project'] == '1') | (METADATA['project'] == '2')].index)

################################
# default rule

rule all:
    input:
        expand('out/unaligned/{sample}.bam', sample = ALL_SAMPLES)

################################
# ensure that fastq files are gzipped

# 50 min per sample
rule gzip_fastq:
    input:
        'data/{file}.fastq'
    output:
        'data/{file}.fastq.gz'
    shell:
        'gzip {input}'

################################
# pipeline

# create unaligned bam file for each sample from the fastq.gz files
# allow ~1 hour per sample for largest samples (15GB of fastq files), 1 core, 3GB memory
rule unaligned_bamfile:
    input:
        read1 = 'data/{sample}_R1.fastq.gz',
        read2 = 'data/{sample}_R2.fastq.gz'
    output:
        'out/unaligned/{sample}.bam'
    conda:
        'envs/picard.yaml'
    #params:
    #    parameter1 = lambda wildcards: SAMPLES.loc[wildcards.sample,'parameter1']
    #threads: 4
    log:
        'out/unaligned/{sample}.bam.log'
    shell:
        'picard FastqToSam F1={input.read1} F2={input.read2} O={output} SM={wildcards.sample} 2>{log}'
