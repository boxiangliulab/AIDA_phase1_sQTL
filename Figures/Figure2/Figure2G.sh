#!/bin/bash

## Run by
##bash 01_run_ggsashimi.sh "chr12:6872841-6873403" "SPSB2" "CD8+_T_GZMB+" 1000
l=$1
g=$2
ct=$3
ind_ethnicity="/ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_ethnicity/list.sample.ethnicity.SG.txt" 
sample_reads=$4



for ethnicity in `cut -f2 $ind_ethnicity | cut -f2 | sort -u`
do
    mkdir -p /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/all_individuals

    # for each individual 
    for ind in `grep $ethnicity $ind_ethnicity | cut -f1`
    do
	samtools view -@ 16 /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/00_raw_bams/$ind.$ct.bam $l \
		 -o /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/all_individuals/$ind.$ct.$ethnicity.$g.bam
    done
    samtools merge -f -@ 16 /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/all_individuals/*$ct.$ethnicity.$g.bam \
	     -o /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.bam
    samtools index -@ 16 /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.bam

    # sample 1000 reads
    reads=$sample_reads
    bam=/ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.bam
    fraction=$(samtools idxstats $bam | cut -f3 | awk -v ct=$reads 'BEGIN {total=0} {total += $1} END {print ct/total}')
    echo $ct
    echo $fraction
    if (( $(echo "$fraction < 1" | bc -l) ))
    then
	echo "samtools view -s ${fraction} /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.bam"
	samtools view -s ${fraction} --subsample-seed 1 /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.bam \
		 -o /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.sampled$sample_reads.bam
    else
	cp /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.bam /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.sampled$sample_reads.bam
    fi
    samtools index -@ 16 /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.sampled$sample_reads.bam

    echo -e "$ethnicity\t/ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/01_gene_ethnicity/$g/mergeAll.$ct.$ethnicity.$g.sampled$sample_reads.bam" >> input_bams.$g.tsv.tmp
done
mv input_bams.$g.tsv.tmp input_bams.$g.tsv


mkdir -p /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/ggsashimi_ethnicity/


~/bin/ggsashimi.py --bam input_bams.$g.tsv \
		   --coordinates $l \
		   --out-prefix /ebs2/aida-sqtl/FreezeV1_0/22_sashimi_plot/ggsashimi_ethnicity/sashimi.$ct.$g.cov50 \
		   --gtf /ebs1/shared/data/reference/hg38/gencode.v32.primary_assembly.annotation.gtf \
		   --palette palette.country.txt -M 50 -C 1 --ann-height 2 \
		   --fix-y-scale --width 8 --height 2 --base-size 18

