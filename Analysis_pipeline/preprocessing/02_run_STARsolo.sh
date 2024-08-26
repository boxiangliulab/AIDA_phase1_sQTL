#!/bin/bash

fullname=$1
abbrname=$2

# run STARsolo
genomeDir=/ebs1/shared/data/reference/STAR_index/hg38_gencode_v32
readFilesIn1=/ebs1/users/tianchi/project/proj_aida_sqtl/data/raw/AIDA_Phase1/$fullname.R1.fastq.gz
readFilesIn2=/ebs1/users/tianchi/project/proj_aida_sqtl/data/raw/AIDA_Phase1/$fullname.R2.fastq.gz
soloCBwhitelist=/ebs1/shared/data/reference/737K-august-2016.txt 
gtfFile=/ebs1/shared/data/reference/hg38/gencode.v32.primary_assembly.annotation.gtf 
vcfFile=/s3/aida-sqtl/FreezeV1_0/01_STARsolo/$abbrname/chrAll.dose.$abbrname.vcf
outDir=/ebs2/aida-sqtl/FreezeV1_0/01_STARsolo/$abbrname

#tmpDir=/ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/01_tmpDir_STAR
#mkdir -p $tmpDir
#rm -r $tmpDir/$fullname

mkdir -p $outDir/log_files

# if you specified --outTmpDir, and this directory exists - please remove it before running STAR
# XS: strand (inferred by STAR based on intron motif), vG: genomic coordinate of the variant, vA: which allele is detected (3 - no match to genotype), vW: WASP filtering 
STAR \
    --runThreadN 8 \
    --genomeDir $genomeDir \
    --twopassMode Basic \
    --soloStrand Forward \
    --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
    --soloBarcodeMate 1 --clip5pNbases 39 0 \
    --soloCBwhitelist $soloCBwhitelist --soloCBmatchWLtype 1MM --soloUMIdedup 1MM_Directional_UMItools \
    --readFilesIn $readFilesIn1 $readFilesIn2 --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM XS vG vA vW \
    --waspOutputMode SAMtag --varVCFfile $vcfFile \
    --sjdbGTFfile $gtfFile \
    --soloOutFileNames solo_$fullname \
    --outFileNamePrefix $outDir/${fullname}_ 1>$outDir/log_files/$fullname.out 2>$outDir/log_files/$fullname.err

#--outTmpDir $tmpDir/$fullname

