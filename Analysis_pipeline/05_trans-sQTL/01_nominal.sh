#!/bin/bash

mkdir -p /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/01_adjust/log/
cd /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/01_adjust/

#/ebs1/zhangyuntian/software/QTLtools/QTLtools trans --vcf genotypes.chr22.vcf.gz --bed genes.simulated.chr22.bed.gz --nominal --normal --out trans.nominal

for celltype in `cat /ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_cell_type/cell_type_use.txt`
do
    for i in `seq 1 22`
    do
	nohup /ebs1/shared/software/qtltools/bin/QTLtools trans \
	      --vcf /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/filter_VCF/chr1_22.sampled.filter.$celltype.dose.map.gt.0.9.vcf.gz \
	      --bed /ebs1/zhangyuntian/project/aida/result/contain_pc_lncRNA/${celltype}${i}_phenotype.txt.gz \
	      --cov /ebs1/zhangyuntian/project/aida/PC_file/new3_21_PC_freeeze/${celltype}_PC.txt \
	      --out /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/01_adjust/${celltype}_pchr${i}_trans.adjust \
	      --adjust /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/permutation/${celltype}_trans.sample.best.txt.gz \
	      --threshold 0.1 --normal --window 5000000 --bin 1000 > /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/01_adjust/log/${celltype}_pchr${i}.log 2>&1 &
	
	
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
    done
    echo "$celltype adjust finished"
done

