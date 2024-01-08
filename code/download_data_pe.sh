#!/bin/bash

# sratoolkit
export PATH=$PATH:$3/bin

# paired-end
fasterq-dump --threads $4 --split-files -O $2 $1

# zip files
gzip $2/$1_1.fastq
gzip $2/$1_2.fastq