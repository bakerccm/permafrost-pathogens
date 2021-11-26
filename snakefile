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
        ['out/fastqc/' + sample for sample in ALL_SAMPLES]

################################

# generate fastqc quality reports for the fastq.gz files in each sample
rule run_fastqc:
    input:
        # list of fastq.gz files for each sample which gets supplied to fastqc command separated by spaces
        read1 = lambda wildcards: METADATA.loc[wildcards.sample,'read1'],
        read2 = lambda wildcards: METADATA.loc[wildcards.sample,'read2']
    output:
        directory('out/fastqc/{sample}')
    conda:
        'envs/fastqc.yaml'
    shell:
        '''
        # remove output directory if it exists; then create new empty directory
            if [ -d {output} ]; then rm -rf {output}; fi
                mkdir {output}
        # run fastqc on each fastq.gz file in the sample
            fastqc -o {output} {input}
        '''

################################

# use multiQC to summarize fastqc results
rule multiQC:
    input:
        fastqc = expand('out/fastqc/{sample}', sample = ALL_SAMPLES) # .zip files from fastqc
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
