### First round: permutation pass

dir=("CD14+_Monocyte" "CD4+_T_em" "CD8+_T_GZMK+" "Treg" "naive_B" "CD16+_Monocyte" "CD4+_T" "CD8+_T_naive"
            "atypical_B" "pDC" "CD16+_NK" "CD4+_T_naive" "IGHMhi_memory_B" "cDC2" "CD4+_T_cm" "CD56+_NK" "IGHMlo_memory_B"
            "gdT_GZMK+" "CD4+_T_cyt" "CD8+_T_GZMB+" "MAIT" "gdT_GZMK-")

##permutation pass
for celltype in ${dir[@]};do
        for i in `seq 1 22`;do
                ./../software/QTLtools/QTLtools cis --vcf ../project/aida/genotype/replication/filter.chr${i}.dose.vcf.gz \
                        --bed ../project/aida/result/contain_pc_lncRNA/${celltype}${i}_phenotype.txt.gz \
                        --cov ../project/aida/PC_file/new3_21_PC_freeeze/${celltype}_PC.txt --permute 1000 --normal --grp-best \
                        --out ../sQTL_result/permute_test_PC/${celltype}_new_${i}_permutation.txt
        done
        cat ../sQTL_result/permute_test_PC/${celltype}_new_*_permutation.txt | gzip -c > ../sQTL_result/permute_test_PC/${celltype}_new_permutations_full.txt.gz
        echo "$celltype permute finished"
done


##FDR correction
for cell in ${dir[@]};do
        Rscript /ebs1/zhangyuntian/software/QTLtools/script/runFDR_cis.R ../permute/${cell}_new_permutations_full.txt.gz  0.05 \
                ../permute/${cell}_new_permutations_full.txt
done

##conditional pass
for cell in ${dir[@]};do
        for i in $(seq 1 22);
        do
        {     ./../../../software/QTLtools/QTLtools cis --vcf ../../../project/aida/genotype/replication/filter.chr${i}.dose.vcf.gz \
                --bed ../../../project/aida/result/contain_pc_lncRNA/${cell}${i}_phenotype.txt.gz \
                --cov ../../../project/aida/PC_file/new3_21_PC_freeeze/${cell}_PC.txt --normal --grp-best \
                --mapping ../permute/${cell}_new_permutations_full.txt.thresholds.txt\
                --out conditions_${cell}_${i}.txt;} &
        done
        wait
        echo "$cell conditional pass"
done


##top variant for each signal

for cell in ${dir[@]};do
        for i in $(seq 1 22);do
                cat conditions_${cell}_${i}.txt | awk '{if($21==1) print $0}' > top_variants_${cell}_${i}.txt
        done
done
