#!/usr/bin/env Rscript

# example R code

#library("here")
#library("tidyverse")
#library("pryr")

# get file names from command line arguments

    args = commandArgs(trailingOnly=TRUE)

    input.filename <- args[1]
    output.filename <- args[2]
    sessioninfo.filename <- args[3]

# code

    # ...
    # ...

# sessionInfo
writeLines(capture.output(sessionInfo()), sessioninfo.filename)
