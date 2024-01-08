#!/bin/bash

conda activate fastp

cd $5

if [ ! -d QCreport ]; then
  mkdir QCreport
fi

##trim adaptor
fastp -i $3/$1 \
	  -o $4/$2_R1_trim.fq.gz \
	  -h $4/QCreport/$2_fastp.html -j $4/QCreport/$2_fastp.json

conda deactivate