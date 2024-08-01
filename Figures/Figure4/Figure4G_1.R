#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Run PCA and Slingshot on B cells 
# author: Yihan Tong
# date: 2023-01-09
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# conda activate r-4.2.3

#' this script is for LMM test of 693 lead sqtl found in B cells: pheno ~ geno + age + sex + 
#' inpput:
#' 1. genotype data: chr*.GT.FORMAT.gz and its variants data: chr1.variants.txt
#' 2. phenotype: 00.index/04.new_phenotype_after_joint.tsv (this is filterd from all sample file according to the sqtl whose phenotyoe is in it)
#' 3. all PCs file
#' 4. sample file with quantile info

library(getopt)
command=matrix(c("file","i",1,"character",
                 "variant","v",1,"character",
                 "chr","c",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if( !is.null(args$help) || is.null(args$file) || is.null(args$variant)|| is.null(args$chr))
{ cat(paste(getopt(command,usage=T),"\n"))
  q() }
geno=args$file #read genotype for each chr#
variant=args$variant
chr=args$chr
chr=as.character(chr)

library(lme4)
library(dplyr)
library(tidyverse)
library(stringr)
#####file and covariants preparation####
# lead sQTL #

sqtl<-read.table("05.intersect_sqtl.only_sqtl.tsv",header=F)
names(sqtl)<-"sqtl"
sqtl$geno<-str_split_fixed(sqtl$sqtl,"_",3)[,3]
sqtl$pheno<-paste0(str_split_fixed(sqtl$sqtl,"_",3)[,1],"_",str_split_fixed(sqtl$sqtl,"_",3)[,2])
sqtl$chr=str_split_fixed(sqtl$sqtl,":",3)[,1]
sqtl<-sqtl[grep(paste0("^", chr, "$"),sqtl$chr),]
# sample with quantile
sample<-read.table("sample_quantile.txt",header = F)
sample$sampleid<-str_split_fixed(sample$V1,"\\.",2)[,1] 
sample$sampleid<-factor(sample$sampleid,levels = unique(sample$sampleid))
sample$pt2<-as.numeric(str_split_fixed(sample$V1,".Q",2)[,2])
sample$pt<-sample$pt2-1
print("01. Finish preparing pseudotime file ^_^")

# genotype file #

genotype<-read.table(geno,sep="\t",header=T)
variant<-read.table(variant,header = F)
colnames(variant)<-"ID"
genotype$ID<-variant
genotype<-genotype[which(genotype$ID$ID %in% sqtl$geno),c("CHROM","POS","ID",sample$V1)] #make the sample order the sample
#geno<-geno[,c(1:2,2621,3:2620)]
##change the type of GT
genotype<-as.data.frame(apply(genotype,2,function(x){gsub(pattern = "0\\|0", replacement = 0,x)}))
genotype<-as.data.frame(apply(genotype,2,function(x){gsub(pattern = "0\\|1", replacement = 1,x)}))
genotype<-as.data.frame(apply(genotype,2,function(x){gsub(pattern = "1\\|0", replacement = 1,x)}))
genotype<-as.data.frame(apply(genotype,2,function(x){gsub(pattern = "1\\|1", replacement = 2,x)}))
print("02. Finish preparing genotype file ^_^")

# phenotype #

phenotype<-read.table("04.new_phenotype_after_joint.tsv",sep="\t",header = T)
phenotype<-phenotype[,c("Chr","start","end","ID","strand","clu",sample$V1)] #make the sample order the sample
print("03. Finish preparing phenotype file ^_^")
# all PCs #
all_pc<-read.table("merge_joint_all_PCs.tsv",sep="\t",header=T)
all_pc<-all_pc[,c("id",sample$V1)] #make the sample order the sample
print("04. Finish preparing PC file ^_^")

#################################################################################################################################################################################


##### Run LMM models #####

# Use lmer to fit the model: intron ratio ~ genotype + pseudotime + genotype*pseudotime + PC

pvalue <- data.frame()
anova_df<-data.frame()
eff<-data.frame()

for (i in 1:nrow(sqtl)) {
  data<-data.frame(sampleid=sample$sampleid,
                   phe=as.numeric(phenotype[which(phenotype$clu==sqtl$pheno[i]),sample$V1]),
                   gen=as.numeric(genotype[which(genotype$ID==sqtl$geno[i]),sample$V1]),
                   sex=as.numeric(all_pc[15,sample$V1]),
                   age=as.numeric(all_pc[14,sample$V1]),
                   # 8 phenotype PC
                   FC1=as.numeric(all_pc[1,sample$V1]),
                   FC2=as.numeric(all_pc[2,sample$V1]),
                   FC3=as.numeric(all_pc[3,sample$V1]),
                   FC4=as.numeric(all_pc[4,sample$V1]),
                   FC5=as.numeric(all_pc[5,sample$V1]),
                   FC6=as.numeric(all_pc[6,sample$V1]),
                   FC7=as.numeric(all_pc[7,sample$V1]),
                   FC8=as.numeric(all_pc[8,sample$V1]),
                   # 5 genotype PC
                   PC1=as.numeric(all_pc[9,sample$V1]),
                   PC2=as.numeric(all_pc[10,sample$V1]),
                   PC3=as.numeric(all_pc[11,sample$V1]),
                   PC4=as.numeric(all_pc[12,sample$V1]),
                   PC5=as.numeric(all_pc[13,sample$V1]),
                   pt=as.numeric(sample$pt))
  data$pt_2<-data$pt*data$pt
  data$geno_pt<-data$gen*data$pt
  data$geno_pt_2<-data$gen*data$pt_2
  # full model #
  fm0 <-lmer(phe ~ gen + pt + geno_pt + FC1 + FC2 + FC3 + FC4 + FC5 + FC6 +FC7 + FC8 + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex + (1 | sampleid),data=data,REML = FALSE)
  coefs <- data.frame(coef(summary(fm0)))
  coefs$p.value <- 2 * pnorm(abs(coefs$t.value), lower.tail = F) # p-value
  #coefs$p.value <- 2 * (1 - pnorm(abs(coefs$t.value)))
  if (c("geno_pt") %in% rownames(coefs)) {
    pvalue_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],p_value=coefs$p.value[which(rownames(coefs)=="geno_pt")])
    eff_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],
                        beta_geno_estimate=as.character(coefs$Estimate[which(rownames(coefs)=="gen")]),
                        beta_inte_estimate=as.character(coefs$Estimate[which(rownames(coefs)=="geno_pt")]))
    eff<-rbind(eff_tmp,eff)
  } else {
    pvalue_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],p_value="none")
  } ### whether "geno_pseudo is dropped"
  rownames(pvalue_tmp)<-NULL
  pvalue<-rbind(pvalue,pvalue_tmp)
  
  # reduced model #
  fm1 <-lmer(phe ~ gen + pt + FC1 + FC2 + FC3 + FC4 + FC5 + FC6 +FC7 + FC8 + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex + (1 | sampleid),data=data,REML = FALSE)
  list<-list(anova = anova(fm1, fm0), fit = fm1)
  anova_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],Pr=(list$anova$`Pr(>Chisq)`[2]))
  anova_df<-rbind(anova_df,anova_tmp)
}
print("05. Finish building Linear mixed models ^_^")
#################################################################################################################################################################################
##### Run quadratic models #####
##use lmer to fit the model: intron ratio ~ genotype + pseudotime + genotype*pseudotime + pseudotime * pseudotime + pseudotime * pseudotime * genotype  + PC
pvalue_qua<-data.frame()
anova_df_qua<-data.frame()
eff_qua<-data.frame()

for (i in 1:nrow(sqtl)) {
	 data<-data.frame(sampleid=sample$sampleid,
                   phe=as.numeric(phenotype[which(phenotype$clu==sqtl$pheno[i]),sample$V1]),
                   gen=as.numeric(genotype[which(genotype$ID==sqtl$geno[i]),sample$V1]),
                   sex=as.numeric(all_pc[15,sample$V1]),
                   age=as.numeric(all_pc[14,sample$V1]),
                   # 8 phenotype PC
                   FC1=as.numeric(all_pc[1,sample$V1]),
                   FC2=as.numeric(all_pc[2,sample$V1]),
                   FC3=as.numeric(all_pc[3,sample$V1]),
                   FC4=as.numeric(all_pc[4,sample$V1]),
                   FC5=as.numeric(all_pc[5,sample$V1]),
                   FC6=as.numeric(all_pc[6,sample$V1]),
                   FC7=as.numeric(all_pc[7,sample$V1]),
                   FC8=as.numeric(all_pc[8,sample$V1]),
                   # 5 genotype PC
                   PC1=as.numeric(all_pc[9,sample$V1]),
                   PC2=as.numeric(all_pc[10,sample$V1]),
                   PC3=as.numeric(all_pc[11,sample$V1]),
                   PC4=as.numeric(all_pc[12,sample$V1]),
                   PC5=as.numeric(all_pc[13,sample$V1]),
                   pt=as.numeric(sample$pt))
  data$pt_2<-data$pt*data$pt
  data$geno_pt<-data$gen*data$pt
  data$geno_pt_2<-data$gen*data$pt_2
  # full model #
  fm0 <-lmer(phe ~ gen + pt + geno_pt + pt_2 + geno_pt_2 + FC1 + FC2 + FC3 + FC4 + FC5 + FC6 +FC7 + FC8 + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex + (1 | sampleid),data=data,REML = FALSE)
  coefs <- data.frame(coef(summary(fm0)))
  coefs$p.value <- 2 * pnorm(abs(coefs$t.value), lower.tail = F) # p-value
  #coefs$p.value <- 2 * (1 - pnorm(abs(coefs$t.value)))
  if (c("geno_pt_2") %in% rownames(coefs)) {
    pvalue_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],p_value=coefs$p.value[which(rownames(coefs)=="geno_pt_2")])
    eff_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],
                        beta_geno_estimate=as.character(coefs$Estimate[which(rownames(coefs)=="gen")]),
                        beta_qua_inte_estimate=as.character(coefs$Estimate[which(rownames(coefs)=="geno_pt_2")]))
    eff_qua<-rbind(eff_tmp,eff_qua)
  } else {
    pvalue_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],p_value="none quadratic")
  } ### whether "geno_pseudo2 is dropped"
  rownames(pvalue_tmp)<-NULL
  pvalue_qua<-rbind(pvalue_qua,pvalue_tmp)
  
  # reduced model #
  fm1 <-lmer(phe ~ gen + pt + pt_2 + geno_pt + FC1 + FC2 + FC3 + FC4 + FC5 + FC6 +FC7 + FC8 + PC1 + PC2 + PC3 + PC4 + PC5 + age + sex + (1 | sampleid),data=data,REML = FALSE)
  list<-list(anova = anova(fm1, fm0), fit = fm1)
  anova_tmp<-data.frame(phenotye=sqtl$pheno[i],variants=sqtl$geno[i],Pr=(list$anova$`Pr(>Chisq)`[2]))
  anova_df_qua<-rbind(anova_df_qua,anova_tmp)
}
print("06. Finish building quadratic models ^_^")

###write lmm results



