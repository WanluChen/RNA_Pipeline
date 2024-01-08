#!/bin/bash

conda activate salmon

salmon quant -i $3 -l A \
         -r $4/$1 \
         -p $6 --validateMappings --gcBias \
         -o $5/$2_quant