#!/bin/bash

mkdir -p /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/log/
cd /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/

#/ebs1/zhangyuntian/software/QTLtools/QTLtools trans --vcf genotypes.chr22.vcf.gz --bed genes.simulated.chr22.bed.gz --nominal --normal --out trans.nominal

for celltype in `cat /ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_cell_type/cell_type_use.txt`
do
    for i in `seq 1 22`
    do
	nohup Rscript ~/project/proj_aida_sqtl/scripts/batch_v1_0/trans_sQTL/runFDR.mappability.R \
	      /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/01_adjust/${celltype}_pchr${i}_trans.adjust.best.txt.gz \
	      0.05 \
	      /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/${celltype}_pchr${i}_trans.adjust.best.genes.txt > /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/log/${celltype}_pchr${i}.log 2>&1 &
	      
	      
	# avoid submit too quickly (pgrep may delay)
	sleep 0.3
	
	# check for number of running jobs constantly, keep running 6 jobs
	while true
	do
            numJobs=`pgrep -cx R -U tianchi`
	    # echo -e "checking # jobs: $numJobs"
            if [ $numJobs -lt 24 ]
            then
		break
            fi
            
            sleep 0.3
	done
    done
    
    echo "$celltype finished"
done

# mappability controled
for celltype in `cat /ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_cell_type/cell_type_use.txt`
do
    echo -e "clu\tpair\tadj_p\tnom_p\tvariant\tgene\tk\tp_min\tp_gene\tqval\tmappability" > /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/all_chr1_22_${celltype}_trans.adjust.best.genes.txt.filtered.txt
    for i in `seq 1 22`
    do
	awk '$11=="PASS"' /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/${celltype}_pchr${i}_trans.adjust.best.genes.txt.filtered.txt >> /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/all_chr1_22_${celltype}_trans.adjust.best.genes.txt.filtered.txt
    done
done

# non controled
for celltype in `cat /ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/summary_cell_type/cell_type_use.txt`
do
    echo -e "clu\tpair\tadj_p\tnom_p\tvariant\tgene\tk\tp_min\tp_gene\tqval" > /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/all_chr1_22_${celltype}_trans.adjust.best.genes.txt
    for i in `seq 1 22`
    do
	cat /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/${celltype}_pchr${i}_trans.adjust.best.genes.txt.filtered.txt >> /ebs1/users/tianchi/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/02_significant_mappability/all_chr1_22_${celltype}_trans.adjust.best.genes.txt
    done
done

