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
        'out/raw/multiqc_report.html', 'out/bbduk/multiqc_report.html', 'out/bbduk_noPhiX/multiqc_report.html', 'out/bbduk_noPhiX_fastuniq/multiqc_report.html' # megahit

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
        read1 = temp('out/bbduk/{sample}_R1.fastq.gz'),
        read2 = temp('out/bbduk/{sample}_R2.fastq.gz'),
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
    threads: 1
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
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################

rule bbduk_noPhiX:
    input:
        read1 = 'out/bbduk/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk/{sample}_R2.fastq.gz'
    output:
        read1_unmatched = temp('out/bbduk_noPhiX/{sample}_unmatched_R1.fastq.gz'), # unmatched reads are not PhiX
        read2_unmatched = temp('out/bbduk_noPhiX/{sample}_unmatched_R2.fastq.gz'),
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
    threads: 1
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
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# deduplicate using FastUniq

# commands for testing
# snakemake --use-conda --rerun-incomplete -j 16 out/bbduk_noPhiX_fastuniq/35m-t0-R2_{R1,R2}.fastq.gz
# snakemake --use-conda --rerun-incomplete -j 56 -np fastuniq_all

rule fastuniq_all:
    input:
        expand('out/bbduk_noPhiX_fastuniq/{sample}_{read}.fastq.gz', sample = ALL_SAMPLES, read = {"R1","R2"})

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
    threads: 28
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
        'multiqc --interactive -o {params.outputdir} {params.inputdir} 2>{log}'

################################
# megahit assembly

# co-assembly of all temperatures and replicates together for each location
# rule megahit:
#     input:
#         read1 = expand("out/bbduk/35m-{temperature}-{replicate}_R1.fastq.gz", temperature = ["t0", "t2"], replicate = ["R1", "R2", "R3"]),
#         read2 = expand("out/bbduk/35m-{temperature}-{replicate}_R2.fastq.gz", temperature = ["t0", "t2"], replicate = ["R1", "R2", "R3"])
#     output:
#         directory("out/megahit")
#     threads: 56
#     conda:
#         'envs/megahit.yaml'
#     shell:
#         '''
#         # get file names and save as array
#             read1_files_space=({input.read1})
#             read2_files_space=({input.read2})
#
#         # use sed to replace spaces with commas (assumes no spaces in file names)
#             read1_files_comma=`echo ${{read1_files_space[@]}} | sed 's/ /,/g'`
#             read2_files_comma=`echo ${{read2_files_space[@]}} | sed 's/ /,/g'`
#
#         # uses printf to join array with commas (works even if file names contain spaces, but these will be passed through and will mess up snakemake unless file names are quoted)
#             # printf -v joined_read1 '%s,' "${{read1_files_space[@]}}"
#             # printf -v joined_read2 '%s,' "${{read2_files_space[@]}}"
#             # read1_files_comma=`echo "${{joined_read1%,}}"`
#             # read2_files_comma=`echo "${{joined_read2%,}}"`
#
#         megahit -1 ${{read1_files_comma}} -2 ${{read2_files_comma}} -t {threads} -o {output}
#         '''

rule megahit:
    input:
        read1 = 'out/bbduk_noPhiX_fastuniq/{sample}_R1.fastq.gz',
        read2 = 'out/bbduk_noPhiX_fastuniq/{sample}_R2.fastq.gz'
    output:
        "out/megahit/{sample}/final.contigs.fa"
    params:
        output_dir = "out/megahit/{sample}"
    threads: 56
    conda:
        'envs/megahit.yaml'
    shell:
        '''
        # remove output directory (megahit will fail if already exists)
            rm -rf {params.output_dir}

        megahit -1 {input.read1} -2 {input.read2} -t {threads} -o {params.output_dir}
        '''

# co-assembles 35m samples, each replicate separately (excluding failed t1 samples)
rule megahit_coassembly_35m_R:
    input:
        read1 = expand('out/bbduk_noPhiX_fastuniq/35m-{temperature}-{{replicate}}_R1.fastq.gz', temperature = ["t0", "t2"]),
        read2 = expand('out/bbduk_noPhiX_fastuniq/35m-{temperature}-{{replicate}}_R2.fastq.gz', temperature = ["t0", "t2"])
    output:
        "out/megahit/35m-{replicate}/final.contigs.fa"
    params:
        output_dir = "out/megahit/35m-{replicate}"
    threads: 56
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

        # remove output directory (megahit will fail if already exists)
            rm -rf {params.output_dir}

        megahit -1 ${{read1_files_comma}} -2 ${{read2_files_comma}} -t {threads} -o {params.output_dir}
        '''

# co-assembles all 35m samples (excluding failed t1 samples)
rule megahit_coassembly_35m:
    input:
        read1 = expand('out/bbduk_noPhiX_fastuniq/35m-{temperature}-{replicate}_R1.fastq.gz', temperature = ["t0", "t2"], replicate = ["R1", "R2", "R3"]),
        read2 = expand('out/bbduk_noPhiX_fastuniq/35m-{temperature}-{replicate}_R2.fastq.gz', temperature = ["t0", "t2"], replicate = ["R1", "R2", "R3"])
    output:
        "out/megahit/35m/final.contigs.fa"
    params:
        output_dir = "out/megahit/35m"
    threads: 56
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

        # remove output directory (megahit will fail if already exists)
            rm -rf {params.output_dir}

        megahit -1 ${{read1_files_comma}} -2 ${{read2_files_comma}} -t {threads} -o {params.output_dir}
        '''

rule megahit_metaquast:
    input:
        "out/megahit/{sample}/final.contigs.fa"
    output:
        directory("out/metaquast/{sample}")
    threads: 4
    conda:
        'envs/quast.yaml'
    shell:
        # use --gene-finding (of -f) to find genes using MetaGeneMark
        # use --max-ref-number 0 to skip searching against SILVA and downloading refs
        'metaquast.py -o {output} -t {threads} --max-ref-number 0 {input}'

################################
