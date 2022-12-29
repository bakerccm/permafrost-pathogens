#!/bin/bash
# Runs in bash on Mac OS X 13.0.1

# Reformat metadata from NCBI

# File names
INPUT="SraRunInfo.csv"
TEMP="samples.temp"
OUTPUT="samples.tsv"

# Samples to exclude from co-assembly
EXCLUDE="35m-t1-R1,35m-t1-R2,35m-t1-R3,83m-t1-R1,83m-t1-R2,83m-t1-R3,CTRL-R1,CTRL-R2,CTRL-R3"

awk -v exclude=${EXCLUDE} \
'BEGIN{
    FS = ","
    OFS="\t"
    # get list of samples to exclude from co-assembly
        split(exclude, exclude_list, ",")
        for (i in exclude_list) exclude_array[exclude_list[i]]
}
NR==1 {
    print "sample", "library_name", "read1", "read2", "location", "time", "temperature_degC", "replicate", "co_assembly"
}
NR!=1 {
    ## all samples ##
    sample = $12
    gsub(/\./, "-", sample)
    library_name = $1
    read1 = $1".sra_1.fastq.gz"
    read2 = $1".sra_2.fastq.gz"
    if (sample ~ /CTRL/)
    {
        ## control samples ##
        location = "NA"
        time = "NA"
        temperature_degC = "NA"
        replicate = "NA"
        co_assembly = "NA"
    }
    else
    {
        ## biological samples ##
        split(sample, sample_characteristics, "-")
        location = sample_characteristics[1]
        time = sample_characteristics[2]
        if (time == "t0")
        {
            temperature_degC = 0
        }
        else if (time == "t1")
        {
            temperature_degC = 3
        }
        else if (time == "t2")
        {
            temperature_degC = 6
        }
        replicate = sample_characteristics[3]
        gsub(/R/, "", replicate)
        if (sample in exclude_array) {co_assembly = "NA"} else {co_assembly = sample_characteristics[1]}
    }
    print sample, library_name, read1, read2, location, time, temperature_degC, replicate, co_assembly
}' ${INPUT} >${TEMP}

# Sort file by sample name
    # Leave header at top
        head -n 1 ${TEMP} >${OUTPUT}
    # Sort remainder of file by sample name
        tail -n +2 ${TEMP} | sort -k 1 >>${OUTPUT}

# Clean up
rm -rf ${TEMP}
unset INPUT TEMP OUTPUT
