###This script was used to generate Figure7B and output enrichment of cis-sQTLs among GWAS associations for 20 immune-mediated or non-immune mediated traits.
###01.prepare annot file
dir=("CD14+_Monocyte" "CD16+_Monocyte" "CD16+_NK" "CD4+_T_cm" "CD4+_T_cyt" "CD4+_T_em" "CD4+_T_naive" \
                 "CD56+_NK" "CD8+_T_GZMB+" "CD8+_T_GZMK+" "CD8+_T_naive" "IGHMhi_memory_B" "IGHMlo_memory_B" "MAIT" "Treg" \
                "cDC2" "gdT_GZMK+" "gdT_GZMK-" "naive_B")

for i in ${dir[@]};do
        for j in `seq 1 22`;do
                python ~/software/ldsc/make_annot.py \
                        --gene-set-file /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/AIDA_EAS_ldscores/gene_annot_file/${i}_gene \
                        --gene-coord-file /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/AIDA_EAS_ldscores/gene_annot_file/gene_pos_aida \
                        --windowsize 100000 \
                        --bimfile /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.${j}.bim \
                        --annot-file /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/AIDA_EAS_ldscores/gene_annot_file/test_${j}.annot.gz
        done
done

###02.compute LD score
for i in ${dir[@]};do
        echo $i
        for j in `seq 1 22`
        do
        {
                echo $j
                python ~/software/ldsc/ldsc.py \
                        --l2  \
                        --bfile /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.${j} \
                        --ld-wind-kb 100 \
                        --annot /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/AIDA_EAS_ldscores/gene_annot_file/${i}_${j}.annot.gz \
                        --thin-annot \
                        --out /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/AIDA_EAS_ldscores/gene_annot_file/${i}_${j} \
                        --print-snps /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/weights/1000G_Phase3_EAS_weights_hm3_no_MHC/test_chr${j}.txt
        } &
        done
        wait
done


###03.compute partition heratibility 
trait=("GBMI_Cell_Genomics_2022_Curated_sumstats_Asthma.txt.gz" \
"HGI_Round7_2022_Curated_sumstats_B2_ALL_eas.txt.gz" \
"Ishigaki_NatGenet_2022_Curated_sumstats_RA_EAS.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_BAS.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_BMI.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_EOS.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_GD.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_Hb.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_Height.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_Ht.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_LYM.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MCH.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MCHC.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MCV.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MON.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_NEU.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_PLT.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_RBC.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_WBC.txt.gz" \
"Shirai_ARD_2022_Curated_sumstats_AD.txt.gz" \
"Shirai_ARD_2022_Curated_sumstats_T1D.txt.gz" \
"Wang_NatCommu_2021_Curated_sumstats_SLE.txt.gz")
celltype=("CD14+_Monocyte" "CD16+_Monocyte" "CD16+_NK" "CD4+_T_cm" "CD4+_T_cyt" "CD4+_T_em" "CD4+_T_naive" \
                 "CD56+_NK" "CD8+_T_GZMB+" "CD8+_T_GZMK+" "CD8+_T_naive" "IGHMlo_memory_B" "MAIT" "Treg" \
                "cDC2" "gdT_GZMK+" "gdT_GZMK-")

for i in ${trait[@]};do
        for j in ${celltype[@]};do
        ~/software/ldsc/ldsc.py \
                --h2 /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/${i}.sumstats.gz \
                --ref-ld-chr /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/AIDA_EAS_ldscores/annot_file/${j}_,/ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/EAS_baselineLD/baseline. \
                --overlap-annot \
                --frqfile-chr /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
                --w-ld-chr /ebs1/zhangyuntian/sQTL_result/colocalization/test/GWAS/test_ldsc/data/weights/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC. \
                --out partition_heratibility_1e_4/${i}_${j} \
                --print-coefficients
        done
done

###barplot (R)
library(ggplot2)
newdat<-read.table("C:/users/90410/desktop/aida/supple_file/table/enrichment_1e-4_all_GWAS.txt",sep="\t",header=TRUE,row.names = 1)
rownames(newdat)<-c("Asthma","RA","BAS","BMI","EOS","GD","Hb","Height","Ht","LYM","MCH","MCHC","MCV","MON","NEU","PLT","RBC","WBC","AD","SLE")
for(i in 1:20){
  newdat[i,20]<-max(as.numeric(newdat[i,1:19]))
}
newdat<-newdat[order(newdat$V20),]
##rearrange newdat
celltype_num<-19
newdata<-as.data.frame(matrix(NA,20*celltype_num,3))
for(i in 1:20){
  for(j in 1:celltype_num){
    newdata[celltype_num*i+j-celltype_num,1]<-rownames(newdat)[i]
    newdata[celltype_num*i+j-celltype_num,2]<-celltype_random[j]
    newdata[celltype_num*i+j-celltype_num,3]<-newdat[i,j]
  }
}
color=c("#fee090","#fdae61", "#bf812d",
        "#9e9ac8","#6a51a3","#807dba",
        "#bcbddc","#8c510a","#3690c0","#74a9cf",
        "#a6bddb", "#a1d99b","#74c476",
        "#0570b0","#CCCC4D","#f46d43",
        "#f768a1", "#dd3497","#c7e9c0")
colnames(newdata)<-c("trait","cell_type","enrichment")
newdata$cell_type<-factor(newdata$cell_type,levels=celltype_random)
newdata$trait<-factor(newdata$trait,levels=rownames(newdat))
ggplot(newdata,aes(x=trait,y=enrichment,fill=cell_type))+geom_bar(position="dodge",stat="identity")+
  scale_fill_manual(values=color)+theme_classic()+coord_flip()+guides(fill=FALSE)
