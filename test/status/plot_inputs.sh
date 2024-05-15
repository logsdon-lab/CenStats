#!/bin/bash

for file in test/input/*.out; do
    bname=$(basename $file .fa.out)
    Rscript test/plot_repeatStructure_onlyRM.R ${file} test/input/${bname}.png
done
