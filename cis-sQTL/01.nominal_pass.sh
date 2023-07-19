dir=("CD14+_Monocyte" "CD16+_Monocyte" "CD16+_NK" "CD4+_T_cm" "CD4+_T_cyt" "CD4+_T_em" "CD4+_T_naive" \
                 "CD56+_NK" "CD8+_T_GZMB+" "CD8+_T_GZMK+" "CD8+_T_naive" "IGHMhi_memory_B" "IGHMlo_memory_B" "MAIT" "Treg" \
                "cDC2" "gdT_GZMK+" "gdT_GZMK-" "naive_B")
for i in `seq 1 22`;do
        for celltype in ${dir[@]}
        do
        {
                        /ebs1/zhangyuntian/software/QTLtools/QTLtools cis --vcf /ebs1/zhangyuntian/project/aida/genotype/replication/filter.chr${i}.dose.vcf.gz \
                        --bed /ebs1/zhangyuntian/project/aida/result/contain_pc_lncRNA/${celltype}${i}_phenotype.txt.gz \
                        --cov /ebs1/zhangyuntian/project/aida/PC_file/new3_21_PC_freeeze/${celltype}_PC.txt --nominal 1 --normal \
                        --grp-best --out /ebs1/zhangyuntian/sQTL_result/nominal_1/${celltype}_${i}_nominals_1.txt
        } &
        done
done
