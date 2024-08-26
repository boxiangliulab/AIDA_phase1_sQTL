#!/bin/bash

grep $1 list_sample.txt | cut -f 2 | sort -u | while read sample
do
    echo $sample
    bash 01_run_prepare_VCF.sh $sample
done
