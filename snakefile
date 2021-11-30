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
        'out/raw/multiqc', 'cutadapt_multiQC', 'sickle_multiQC'

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
        'fastqc -o out/raw/fastqc {input}'

# use multiQC to summarize fastqc results
rule raw_multiQC:
    input:
        fastqc = expand('out/raw/fastqc/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
        # star = expand('out/alignment/{sample}_Log.final.out', sample = list(SAMPLES.index)), # .Log.final.out files from STAR alignment
        # featureCounts = 'out/featureCounts/bulkseq_featureCounts.txt.summary', # .summary file from featureCounts
        # htseq_count = expand('out/htseq_count/{sample}.txt', sample = list(SAMPLES.index)) # .txt files from htseq-count
    output:
        directory('out/raw')
    params:
        inputdir = 'out/raw'
    log:
        'out/raw/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc -o {output} {params.inputdir} 2>{log}'

################################
# remove sequencing adaptors

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
        adapter_fwd = config['cutadapt']['adapter_fwd'],
        adapter_rev = config['cutadapt']['adapter_rev'],
        min_length = config['cutadapt']['min_length']
    conda:
        'envs/cutadapt-3.5.yaml'
    threads: 4
    shell:
        '''
        cutadapt -a {params.adapter_fwd} -A {params.adapter_rev} \
        -j {threads} -m {params.min_length} \
        -o {output.read1} -p {output.read2} {input.read1} {input.read2} >{output.qc}
        '''

################################

rule cutadapt_fastqc:
    input:
        'out/cutadapt/{sample}_{read}.fastq.gz'
    output:
        'out/cutadapt/fastqc/{sample}_{read}_fastqc.html',
        'out/cutadapt/fastqc/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/cutadapt/fastqc {input}'

rule cutadapt_multiQC:
    input:
        fastqc = expand('out/cutadapt/fastqc/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
        # star = expand('out/alignment/{sample}_Log.final.out', sample = list(SAMPLES.index)), # .Log.final.out files from STAR alignment
        # featureCounts = 'out/featureCounts/bulkseq_featureCounts.txt.summary', # .summary file from featureCounts
        # htseq_count = expand('out/htseq_count/{sample}.txt', sample = list(SAMPLES.index)) # .txt files from htseq-count
    output:
        directory('out/cutadapt/multiqc')
    params:
        inputdir = 'out/cutadapt/fastqc'
    log:
        'out/cutadapt/multiqc/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc -o {output} {params.inputdir} 2>{log}'

################################
## trim reads using sickle

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
    # log:
    #     'out/sickle/{sample}.log'
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
            -q {params.quality_threshold} -l {params.length_threshold}
        '''

rule sickle_fastqc:
    input:
        'out/sickle/{sample}_{read}.fastq.gz'
    output:
        'out/sickle/fastqc/{sample}_{read}_fastqc.html',
        'out/sickle/fastqc/{sample}_{read}_fastqc.zip'
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc -o out/sickle/fastqc {input}'

rule sickle_multiQC:
    input:
        fastqc = expand('out/sickle/fastqc/{sample}_{read}_fastqc.zip', sample = ALL_SAMPLES, read = {'R1','R2'})
        # star = expand('out/alignment/{sample}_Log.final.out', sample = list(SAMPLES.index)), # .Log.final.out files from STAR alignment
        # featureCounts = 'out/featureCounts/bulkseq_featureCounts.txt.summary', # .summary file from featureCounts
        # htseq_count = expand('out/htseq_count/{sample}.txt', sample = list(SAMPLES.index)) # .txt files from htseq-count
    output:
        directory('out/sickle/multiqc')
    params:
        inputdir = 'out/sickle/fastqc'
    log:
        'out/sickle/multiqc/multiqc_report.log'
    conda:
        'envs/multiqc.yaml'
    shell:
        'multiqc -o {output} {params.inputdir} 2>{log}'

################################
## merge reads using fastq-join in ea-utils
# see https://github.com/ExpressionAnalysis/ea-utils/blob/wiki/FastqJoin.md

rule fastq_join:
    input:
        read1 = 'out/sickle/{sample}_R1.fastq.gz',
        read2 = 'out/sickle/{sample}_R2.fastq.gz'
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
        -o {output.join} {output.unmatched1} {output.unmatched2}
        '''

################################
