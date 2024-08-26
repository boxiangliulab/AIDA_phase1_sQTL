#!/bin/bash

mkdir -p /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/permutation/log
cd /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/permutation

#/ebs1/zhangyuntian/software/QTLtools/QTLtools trans --vcf genotypes.chr22.vcf.gz --bed genes.simulated.chr22.bed.gz --nominal --normal --out trans.nominal

for celltype in `cat /ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_cell_type/cell_type_use.txt`
do
    nohup QTLtools trans \
	  --vcf /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/filter_VCF/chr1_22.sampled.filter.$celltype.dose.map.gt.0.9.vcf.gz \
	  --bed /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/merged_bed/chr1_22.${celltype}_phenotype.txt.gz \
	  --cov /ebs1/zhangyuntian/project/aida/PC_file/new3_21_PC_freeeze/${celltype}_PC.txt \
	  --out /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/permutation/${celltype}_trans.sample \
	  --sample 50000 \
	  --threshold 1e-5 --normal --window 5000000 --bin 1000 > /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/permutation/log/${celltype}.log 2>&1 &


    # avoid submit too quickly (pgrep may delay)
    sleep 1
    
    # check for number of running jobs constantly, keep running 6 jobs
    while true
    do
        numJobs=`pgrep -cx QTLtools -U tianchi`
	# echo -e "checking # jobs: $numJobs"
        if [ $numJobs -lt 32 ]
        then
	    break
        fi
        
        sleep 1
    done
    
    echo "$celltype permutation finished"
done

