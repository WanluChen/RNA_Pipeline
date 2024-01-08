#!/bin/bash

conda activate fastp

cd $5

if [ ! -d QCreport ]; then
  mkdir QCreport
fi

##trim adaptor
fastp -i $4/$1 -I $4/$2 \
	  -o $5/$3_R1_trim.fq.gz -O $5/$3_R2_trim.fq.gz \
	  -h $5/QCreport/$3_fastp.html -j $5/QCreport/$3_fastp.json

conda deactivate