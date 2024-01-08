#!/bin/bash

# sratoolkit
export PATH=$PATH:$3/bin

# single-end
fasterq-dump --threads $4 -O $2 $1

# zip files
gzip $2/$1.fastq