###This script was used to identify sex-biased sQTL from independent cis-sQTLs we found previously.
library(qvalue)
library(data.table)
##prepare genotype file

genocol<-read.table("/ebs1/zhangyuntian/script/sex_biased_sQTL/geno.txt",sep="\t")
celltype = c("CD14+_Monocyte","CD16+_Monocyte","CD16+_NK","CD4+_T_cm","CD4+_T_cyt","CD4+_T_em","CD4+_T_naive",
                  "CD56+_NK","CD8+_T_GZMB+","CD8+_T_GZMK+","CD8+_T_naive","IGHMhi_memory_B","IGHMlo_memory_B","MAIT","Treg",
                  "cDC2","gdT_GZMK+","gdT_GZMK-","naive_B")

for(cell in 1:length(celltype)){
        condition<-as.data.frame(matrix(NA,0,22))
  for(j in 1:22){
          file_name=paste("/ebs1/zhangyuntian/sQTL_result/conditional_pass/new_3_23_condition_result/top_variants_",celltype[cell],"_",j,".txt",sep="")
          if(file.info(file_name)$size != 0){
          condition_sQTL<-read.table(paste("/ebs1/zhangyuntian/sQTL_result/conditional_pass/new_3_23_condition_result/top_variants_",celltype[cell],"_",j,".txt",sep=""),sep=" ")
          condition<-rbind(condition,condition_sQTL)}}
  print("finish read in condition result")
##read in genotype file
  genotype<-as.data.frame(matrix(NA,0,ncol(genocol)))
  for(j in 1:nrow(condition)){
          geno_per_SNP<-read.table(paste("/ebs1/zhangyuntian/script/sex_biased_sQTL/SNP_info/",condition[j,10],".txt",sep=""),sep="\t")
          genotype<-rbind(genotype,geno_per_SNP)
  }
  colnames(genotype)<-genocol[1,]
##intersect
  geno_cell<-genotype
  geno_cell<-geno_cell[,-c(1:9)]
  for(i in 1:nrow(geno_cell)){
    for(j in 1:ncol(geno_cell)){
      if(startsWith(geno_cell[i,j],'0|0')){geno_cell[i,j]<-0}
      if(startsWith(geno_cell[i,j],'0|1')){geno_cell[i,j]<-1}
      if(startsWith(geno_cell[i,j],'1|0')){geno_cell[i,j]<-1}
      if(startsWith(geno_cell[i,j],'1|1')){geno_cell[i,j]<-2}
    }}
  ##prepare sex file & PC file
  PC<-read.table(paste("/ebs1/zhangyuntian/project/aida/PC_file/new3_21_PC_freeeze/",celltype[cell],"_PC.txt",sep=""),sep="\t",header=TRUE,row.names=1)
  ##prepare phenotype file
  col_phe<-read.table(paste("/ebs1/zhangyuntian/script/sex_biased_sQTL/phename/",celltype[cell],"_phename.txt",sep=""),sep=" ")
  newphe<-as.data.frame(matrix(NA,0,ncol(col_phe)))
  for(j in 1:22){
          phenotype<-fread(paste("/ebs1/zhangyuntian/project/aida/result/contain_pc_lncRNA/",celltype[cell],j,"_phenotype.txt.gz",sep=""),sep="\t")
          phenotype<-as.data.frame(phenotype)
          file_name=paste("/ebs1/zhangyuntian/sQTL_result/conditional_pass/new_3_23_condition_result/top_variants_",celltype[cell],"_",j,".txt",sep="")
          if(file.info(file_name)$size != 0){
          condition_sQTL<-read.table(paste("/ebs1/zhangyuntian/sQTL_result/conditional_pass/new_3_23_condition_result/top_variants_",celltype[cell],"_",j,".txt",sep=""),sep=" ")
          for(ro in 1:nrow(condition_sQTL)){
                  num<-which(phenotype[,4]==condition_sQTL[ro,6])
                  newphe<-rbind(newphe,phenotype[num,])
          }
          print(j)}
  }
  phenotype<-newphe
  phenotype<-phenotype[,-c(1:6)]
  inter<-intersect(colnames(geno_cell),colnames(phenotype))
  inter<-intersect(colnames(PC),inter)
  phenotype<-phenotype[,inter]
  geno_cell<-geno_cell[,inter]
  PC<-PC[,inter]


##use lmer to fit the model: intron ratio ~ genotype + sex + genotype | sex + PC

  pvalue<-as.data.frame(matrix(1,nrow(condition),5))

  for(i in 1:nrow(condition)){
    phe<-as.numeric(t(phenotype[i,]))
    geno<-as.numeric(t(geno_cell[i,]))
    pvalue[i,1]<-condition[i,6]
    pvalue[i,2]<-condition[i,10]
    PC1<-as.numeric(t(PC[1,]))
    PC2<-as.numeric(t(PC[2,]))
    PC3<-as.numeric(t(PC[3,]))
    PC4<-as.numeric(t(PC[4,]))
    PC5<-as.numeric(t(PC[5,]))
    PC6<-as.numeric(t(PC[6,]))
    PC7<-as.numeric(t(PC[7,]))
    PC8<-as.numeric(t(PC[8,]))
    PC9<-as.numeric(t(PC[9,]))
    PC10<-as.numeric(t(PC[10,]))
    PC11<-as.numeric(t(PC[11,]))
    PC12<-as.numeric(t(PC[12,]))
    PC13<-as.numeric(t(PC[13,]))
    sex<-as.numeric(t(PC[14,]))
    PC15<-as.numeric(t(PC[15,]))
    geno_sex<-geno*sex
    fm <-lm(phe ~ geno + sex + geno*sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC15)
    pvalue[i,3]<-summary(fm)$coefficients[,'Pr(>|t|)']['geno:sex']
    pvalue[i,4]<-summary(fm)$coefficients[,'Estimate']['geno:sex']
    pvalue[i,5]<-summary(fm)$coefficients[,'Std. Error']['geno:sex']
  }
   write.table(pvalue,paste("pvalue",celltype[cell],".txt",sep=""),sep="\t",quote=FALSE)
   print("pvalue computed!")
   pvalue<-na.omit(pvalue)
   qvalue<-pvalue
   for(i in 1:nrow(qvalue)){
     name<-pvalue[i,1]
     num<-which(pvalue[,1]==name)
     p<-pvalue[num,3]
     qvalue[num,3]<-p.adjust(p,method="bonferroni",n=length(p))}
     print("bonferroni correction finished")
   length(which(qvalue[,3]<0.25))
   storeyqvalue<-as.data.frame(matrix(1,nrow(qvalue),3))
   for(i in 1:nrow(qvalue)){
     storeyqvalue[i,1]<-qvalue[i,1]
     storeyqvalue[i,2]<-qvalue[i,2]
     name<-qvalue[i,1]
     num<-which(qvalue[,1]==name)
     least<-min(qvalue[num,3])
     if(least!=qvalue[i,3])storeyqvalue[i,2]<-NA
     storeyqvalue[i,3]<-least }
   newqvalue<-na.omit(storeyqvalue)
   final_qvalue<-qvalue(newqvalue[,3],fdr.level=0.25,pi0=1)$qvalues
   print("FDR correction finished")
   newqvalue<-cbind(newqvalue,final_qvalue)
   write.table(newqvalue,paste("final_qvalue_",celltype[cell],".txt",sep=""),sep="\t",quote=FALSE)}

###compute interaction effect size & p-value 
library(data.table)
##compute beta & pvalue per sex
beta_pvalue<-as.data.frame(matrix(NA,89,7))
genocol<-read.table("/ebs1/zhangyuntian/script/sex_biased_sQTL/geno.txt",sep="\t")

celltype = c("CD14+_Monocyte","CD16+_Monocyte","CD16+_NK","CD4+_T_cm","CD4+_T_em","CD4+_T_naive",
                  "CD8+_T_GZMB+","CD8+_T_GZMK+","CD8+_T_naive","IGHMhi_memory_B","IGHMlo_memory_B","MAIT",
                  "cDC2","gdT_GZMK-")
check=0
for(i in 1:length(celltype)){
        print(i)
        intron_snp<-read.table(paste("intron_snp_sbQTL",celltype[i],".txt",sep=""),sep="\t")
        ##read in phenotype
        col_phe<-read.table(paste("/ebs1/zhangyuntian/script/sex_biased_sQTL/phename/",celltype[i],"_phename.txt",sep=""),sep=" ")
        newphe<-as.data.frame(matrix(NA,0,ncol(col_phe)))
        for(j in 1:22){
                phenotype<-fread(paste("/ebs1/zhangyuntian/project/aida/result/contain_pc_lncRNA/",celltype[i],j,"_phenotype.txt.gz",sep=""),sep="\t")
                phenotype<-as.data.frame(phenotype)
                for(ro in 1:nrow(intron_snp)){
                        num<-which(phenotype[,4]==intron_snp[ro,2])
                        newphe<-rbind(newphe,phenotype[num,])
                }
        }
        phenotype<-newphe
        PC<-read.table(paste("/ebs1/zhangyuntian/project/aida/PC_file/new3_21_PC_freeeze/",celltype[i],"_PC.txt",sep=""),sep="\t",header=TRUE,row.names=1)
        for(k in 1:nrow(intron_snp)){
                genotype<-read.table(paste("SNP_info/",intron_snp[k,3],".txt",sep=""),sep="\t")
                colnames(genotype)<-genocol[1,]
                num<-which(phenotype[,4]==intron_snp[k,2])
                pheno<-phenotype[num,]
                pheno<-pheno[,-c(1:6)]
                inter<-intersect(colnames(genotype),colnames(pheno))
                inter<-intersect(colnames(PC),inter)
                pheno<-as.data.frame(pheno)
                pheno<-pheno[,inter]
                genotype<-as.data.frame(genotype)
                genotype<-genotype[,inter]
                geno_cell<-genotype
                for(m in 1:nrow(geno_cell)){
                        for(n in 1:ncol(geno_cell)){
                                if(startsWith(geno_cell[m,n],'0|0')){geno_cell[m,n]<-0}
                                if(startsWith(geno_cell[m,n],'0|1')){geno_cell[m,n]<-1}
                                if(startsWith(geno_cell[m,n],'1|0')){geno_cell[m,n]<-1}
                                if(startsWith(geno_cell[m,n],'1|1')){geno_cell[m,n]<-2}}}
                PC<-PC[,inter]
                sex<-PC[14,]
                male_num<-which(sex[1,]==0)
                female_num<-which(sex[1,]==1)
                male_pheno<-pheno[,colnames(sex)[male_num]]
                male_geno<-geno_cell[,colnames(sex)[male_num]]

                female_pheno<-pheno[,colnames(sex)[female_num]]
                female_geno<-geno_cell[,colnames(sex)[female_num]]

                male_fm <-lm(as.numeric(t(male_pheno)) ~ as.numeric(t(male_geno)))
                female_fm<- lm(as.numeric(t(female_pheno)) ~ as.numeric(t(female_geno)))
                check = check + 1
                beta_pvalue[check,1]<-celltype[i]
                beta_pvalue[check,2]<-intron_snp[k,3]
                beta_pvalue[check,3]<-intron_snp[k,2]
                beta_pvalue[check,4]<-summary(male_fm)$coefficients[,'Estimate']['as.numeric(t(male_geno))']
                beta_pvalue[check,5]<-summary(male_fm)$coefficients[,'Pr(>|t|)']['as.numeric(t(male_geno))']
                beta_pvalue[check,6]<-summary(female_fm)$coefficients[,'Estimate']['as.numeric(t(female_geno))']
                beta_pvalue[check,7]<-summary(female_fm)$coefficients[,'Pr(>|t|)']['as.numeric(t(female_geno))']
        }
}

write.table(beta_pvalue,"beta_pvalue_sbQTL.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
