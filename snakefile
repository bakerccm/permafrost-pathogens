# snakefile for permafrost pathogens processing and analysis

# Chris Baker
# https://github.com/bakerccm
# 26 November 2021

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

################################
# default rules

rule all:
    input:
        'out/raw/multiqc_report.html', 'out/bbduk/multiqc_report.html', 'out/bbduk_noPhiX/multiqc_report.html', 'out/fastq-join/multiqc_report.html'

################################
rule all_raw_data_links:
    input:
        expand('data/links/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {'R1','R2'})

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
# see bbduk overview here http://seqanswers.com/forums/showthread.php?t=42776
rule bbduk:
    input:
        read1 = 'data/links/{sample}_R1.fastq.gz',
        read2 = 'data/links/{sample}_R2.fastq.gz'
    output:
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
    shell:
        '''
        bbduk.sh {params.memory} in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} \
        # adapter filtering
            ref={params.ref} ktrim={params.ktrim} k={params.k} mink={params.trim_params} hdist={params.hdist} {params.trim_params} \
        # quality filtering
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
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################

rule bbduk_noPhiX:
    input:
        read1 = 'out/bbduk/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk/{sample}_R2.fastq.gz'
    output:
        read1_unmatched = 'out/bbduk_noPhiX/{sample}_unmatched_R1.fastq.gz', # unmatched reads are not PhiX
        read2_unmatched = 'out/bbduk_noPhiX/{sample}_unmatched_R2.fastq.gz',
        read1_matched = 'out/bbduk_noPhiX/{sample}_matched_R1.fastq.gz', # matched reads are PhiX; consider just discarding these
        read2_matched = 'out/bbduk_noPhiX/{sample}_matched_R2.fastq.gz',
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
    shell:
        '''
        bbduk.sh {params.memory} in1={input.read1} in2={input.read2} out1={output.read1_unmatched} out2={output.read2_unmatched} outm1={output.read1_matched} outm2={output.read2_matched} \
        # contaminant filtering
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
        read1 = expand("out/bbduk/35m-{temperature}-{replicate}_R1.fastq.gz", temperature = ["t0", "t2"], replicate = ["R1", "R2", "R3"]),
        read2 = expand("out/bbduk/35m-{temperature}-{replicate}_R2.fastq.gz", temperature = ["t0", "t2"], replicate = ["R1", "R2", "R3"])
    output:
        directory("out/megahit")
    threads: 56
    conda:
        'envs/megahit.yaml'
    shell:
        '''
        # get file names and save as array
            read1_files_space=({input.read1})
            read2_files_space=({input.read2})
        
        # use sed to replace spaces with commas (assumes no spaces in file names)
            read1_files_comma=`echo ${{read1_files_space[@]}} | sed 's/ /,/g'`
            read2_files_comma=`echo ${{read2_files_space[@]}} | sed 's/ /,/g'`
        
        # uses printf to join array with commas (works even if file names contain spaces, but these will be passed through and will mess up snakemake unless file names are quoted)
            # printf -v joined_read1 '%s,' "${{read1_files_space[@]}}"
            # printf -v joined_read2 '%s,' "${{read2_files_space[@]}}"
            # read1_files_comma=`echo "${{joined_read1%,}}"`
            # read2_files_comma=`echo "${{joined_read2%,}}"`
        
        megahit -1 ${{read1_files_comma}} -2 ${{read2_files_comma}} -t {threads} -o {output}
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
    threads: 32
    conda:
        'envs/spades.yaml'
    shell:
        # note hard codes memory limit of 180Gb here (Frontera nodes have 192Gb so this should run on a single node if it's the only job running)
        # note that this fails OOM if the limit is set to 190Gb
        '''
        metaspades.py -t {threads} -m 180 \
        --pe1-1 {input.pe1_1} --pe1-2 {input.pe1_2} \
        --pe2-1 {input.pe2_1} --pe2-2 {input.pe2_2} \
        -o {output}
        '''

################################
