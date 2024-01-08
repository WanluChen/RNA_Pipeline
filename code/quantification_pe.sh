#!/bin/bash

conda activate salmon

salmon quant -i $4 -l A \
         -1 $5/$1 \
         -2 $5/$2 \
         -p $7 --validateMappings --gcBias \
         -o $6/$3_quant