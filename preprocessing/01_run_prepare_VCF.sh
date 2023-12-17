#!/bin/bash

sample=$1

# prepare VCF
#bcftools concat `ls /ebs1/shared/data/aida/Imputation/PostIMP/chr*.dose.vcf.gz` --threads 4 -o /s3/aida-sqtl/FreezeV1_0/01_STARsolo/chrAll.dose.vcf.gz
mkdir -p /s3/aida-sqtl/FreezeV1_0/01_STARsolo/$sample/
bcftools view --threads 4 /s3/aida-sqtl/FreezeV1_0/01_STARsolo/chrAll.dose.vcf.gz -s $(grep $sample conversion_DCPID_sampleID.csv | cut -d' ' -f1) -o /s3/aida-sqtl/FreezeV1_0/01_STARsolo/$sample/chrAll.dose.$sample.vcf.gz
gunzip /s3/aida-sqtl/FreezeV1_0/01_STARsolo/$sample/chrAll.dose.$sample.vcf.gz
