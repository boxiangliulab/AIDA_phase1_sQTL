#!/bin/bash

mkdir -p heatmap
cd heatmap

samtools index -@ 16 ../mergeAll001.all_celltypes.R1.bam
samtools index -@ 16 ../mergeAll001.all_celltypes.R2.bam

bamCoverage -b ../mergeAll001.all_celltypes.R1.bam -o ../mergeAll001.all_celltypes.R1.bw --normalizeUsing RPKM -bs 5 -p 16 >log.bw1.txt 2>&1
bamCoverage -b ../mergeAll001.all_celltypes.R2.bam -o ../mergeAll001.all_celltypes.R2.bw --normalizeUsing RPKM -bs 5 -p 16 >log.bw2.txt 2>&1


# final version: expressed genes
grep -f ../genelist_expressed.txt /ebs1/shared/data/reference/hg38/gencode.v32.annotation.gtf | less > gencode.v32.annotation.expr.gtf
computeMatrix scale-regions \
              -R gencode.v32.annotation.expr.gtf \
              -S ../mergeAll001.all_celltypes.R1.bw ../mergeAll001.all_celltypes.R2.bw \
              --metagene \
              -o mergeAll001.all_celltypes.mat.final.gz --outFileSortedRegions sortedRegions.bed \
              -b 1000 -a 1000 \
              --sortUsing mean \
              --skipZeros \
              -p 16 >log.heatmap.matrix.final.txt 2>&1

plotHeatmap -m mergeAll001.all_celltypes.mat.final.gz \
            --samplesLabel Read1 Read2 \
            --heatmapWidth 4 --heatmapHeight 9 \
            --colorMap Blues Blues \
            -out mergeAll001.all_celltypes.heatmap.final.pdf >log.heatmap.final.txt 2>&1
