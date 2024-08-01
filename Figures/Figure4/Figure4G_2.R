library(ggplot2)
library(pheatmap)
library(patchwork)
library(tidyverse)
library(stringr)
library(dplyr)
library(pheatmap)

delete.na<-function(df,n=0) {
  df[rowSums(is.na(df)) <= n,]
}
####plot LMMs anova test value
#anova_lmm<-read.csv("/Users/tongyihan/project/AIDA/data/dynamic sQTL/test JP/test_v2/LMMs_anova_test.csv",row.names = 1)
anova_lmm<-read.table("merge.LMM.anova_test_all.tsv",header = T)
anova_qua<-read.table("merge.quadratic.anova_test_all.tsv",header = T)
anova_lmm<-delete.na(anova_lmm)

###plot phenotype value (for those sqtls with p_adi <= 0.05)
#####################################################################################################################################################################################################################
###read files

#phenotype
#phenotype<-read.csv("/Users/tongyihan/project/AIDA/data/dynamic sQTL/test JP/JP_test_perind.counts.gz.qqnorm_chr22",sep="\t")
phenotype<-read.csv("JP_test_perind.counts.gz.qqnorm_all_chr",sep="\t")
phenotype<-read.csv("04.new_phenotype_after_joint.tsv",sep="\t")
phenotype<-phenotype[,-c(1:4)]
phenotype<-phenotype[,-2619]
#phenotype<-dplyr::filter(phenotype,ID %in% intersect(phenotype$ID,nominal_sQTL$V1))
phenotype<-column_to_rownames(phenotype,var="clu")
#phenotype<-phenotype[,-grep("H121",colnames(phenotype))]
print("finish preparing phenotype data ^_^")


lmm_sqtl<-read.table("merge.LMM.anova_test_0.05.tsv",header = T)
lmm_sqtl$model<-"linear model"
lmm_sqtl$sqtl<-paste0(lmm_sqtl$phenotye,"_",lmm_sqtl$variants)
qua_sqtl<-read.table("merge.quadratic.anova_test_0.05.tsv",header = T)
qua_sqtl$model<-"quadratic model"
qua_sqtl$sqtl<-paste0(qua_sqtl$phenotye,"_",qua_sqtl$variants)


#genotype
genotype<-read.table("merge_all_chr_unique_geno.tsv",header = T)
genotype<-genotype[,-c(1:2)]
genotype<-column_to_rownames(genotype,var = "ID")

#####################################################################################################################################################################################################################
intersect<-intersect(lmm_sqtl$sqtl,qua_sqtl$sqtl)
intersect_qtl<-lmm_sqtl[which(lmm_sqtl$sqtl %in% intersect(lmm_sqtl$sqtl,qua_sqtl$sqtl)),]
intersect_qtl$model="Linear & Quadratic"
only_qua<-qua_sqtl[which(qua_sqtl$sqtl%in%setdiff(qua_sqtl$sqtl,lmm_sqtl$sqtl)),]
only_lmm<-lmm_sqtl[which(lmm_sqtl$sqtl%in%setdiff(lmm_sqtl$sqtl,qua_sqtl$sqtl)),]

##effect size pheatmap
###################################the previous one: all linear#################################
lmm_eff<-read.table("merge.LMM.effect_size.tsv",header = T)
lmm_eff$sqtl<-paste0(lmm_eff$phenotye,"_",lmm_eff$variants)
lmm_eff<-lmm_eff[which(lmm_eff$sqtl %in% lmm_sqtl$sqtl),]
lmm_eff$effect_size<-"a"
lmm_heatmap_plot<-data.frame(sqtl=lmm_eff$sqtl,
                             Q1=lmm_eff$beta_geno_estimate+0*lmm_eff$beta_inte_estimate,
                             Q2=lmm_eff$beta_geno_estimate+1*lmm_eff$beta_inte_estimate,
                             Q3=lmm_eff$beta_geno_estimate+2*lmm_eff$beta_inte_estimate,
                             Q4=lmm_eff$beta_geno_estimate+3*lmm_eff$beta_inte_estimate,
                             Q5=lmm_eff$beta_geno_estimate+4*lmm_eff$beta_inte_estimate,
                             Q6=lmm_eff$beta_geno_estimat+5*lmm_eff$beta_inte_estimate)


qua_eff<-read.table("merge.quadratic.effect_size.tsv",header = T)
qua_eff$sqtl<-paste0(qua_eff$phenotye,"_",qua_eff$variants)
qua_eff<-qua_eff[which(qua_eff$sqtl %in% qua_sqtl$sqtl),]
qua_heatmap_plot<-data.frame(sqtl=qua_eff$sqtl,
                             Q1=qua_eff$beta_geno_estimate+0*qua_eff$beta_qua_inte_estimate,
                             Q2=qua_eff$beta_geno_estimate+1*qua_eff$beta_qua_inte_estimate,
                             Q3=qua_eff$beta_geno_estimate+2*qua_eff$beta_qua_inte_estimate,
                             Q4=qua_eff$beta_geno_estimate+3*qua_eff$beta_qua_inte_estimate,
                             Q5=qua_eff$beta_geno_estimate+4*qua_eff$beta_qua_inte_estimate,
                             Q6=qua_eff$beta_geno_estimat+5*qua_eff$beta_qua_inte_estimate)

# the plot mapped together
#merge_heatmap<-rbind(lmm_heatmap_plot,qua_heatmap_plot)
#merge_heatmap<-column_to_rownames(merge_heatmap,var="sqtl")
#merge_heatmap<-merge_heatmap %>% group_by(sqtl) %>% mutate(count=row_number())
#merge_heatmap<-merge_heatmap[,-8]
#merge_heatmap<-column_to_rownames(merge_heatmap,var = "sqtl")
#pheatmap(merge_heatmap,cluster_cols = F,annotation_row = annotation1,clustering_distance_rows = "euclidean")

#for (i in 1:nrow(merge_heatmap)) {
#  if (merge_heatmap$count[i]>1) {
#    merge_heatmap$sqtl[i]=paste0(merge_heatmap$sqtl[i],"_",merge_heatmap$count[i])
#  }
#}

###################################new one recall effect size using QTLtools separately#################################
library(ComplexHeatmap)
library(RColorBrewer)
eff_size<-read.table("merge.beta.tsv",sep="\t",header = T)

# linear heatmap
lmm_heatmap_plot2<-eff_size[which(eff_size$sqtl %in% lmm_sqtl2$sqtl),]
rownames(lmm_heatmap_plot2)<-NULL
lmm_heatmap_plot2<-column_to_rownames(lmm_heatmap_plot2,var = "sqtl")
lmm_heatmap_plot_scale2<-t(scale(t(lmm_heatmap_plot2),scale = F))

lmm_sqtl2<-lmm_sqtl
for (i in which(lmm_sqtl2$sqtl %in% intersect)) {
  lmm_sqtl2$model[i]="Linear & Quadratic"
  lmm_sqtl2$sqtl[i]<-paste0(lmm_sqtl2$sqtl[i],"_2")
}
# the plot that separated
lmm_sqtl2<-lmm_sqtl
for (i in which(lmm_sqtl2$sqtl %in% intersect)) {
  lmm_sqtl2$model[i]="Linear & Quadratic"
}

#plot1-lmm
anno_col<-list(model1=c("Linear & Quadratic"="#e78ac3","linear model"="#a6d854"),
               model2=c("linear model"="#a6d854"))

#annotation for heatmap
annotation_lmm1<-lmm_sqtl2[,c("sqtl","model")]
annotation_lmm2<-lmm_sqtl[,c("sqtl","model")]
annotation_lmm<-merge(annotation_lmm1,annotation_lmm2,by="sqtl")
names(annotation_lmm)<-c("sqtl","model1","model2")
annotation_lmm<-column_to_rownames(annotation_lmm,var="sqtl")
#lmm_heatmap_plot<-column_to_rownames(lmm_heatmap_plot,var="sqtl")

#annotation for complexheatmap
col_anno<-c("Linear & Quadratic"="#66a61e","linear model"="#e6ab02")
annotation_lmm$color<-col_anno[annotation_lmm$model1]
annotation_lmm$sqtl <- rownames(annotation_lmm)
ha = rowAnnotation(bar = annotation_lmm$sqtl,
                   col = list(bar=setNames(as.character(annotation_lmm$color), as.character(annotation_lmm$sqtl))),
                   show_legend=c(bar=FALSE))

pdf("Fig4.LMM.complexheatmap_scale_v2.pdf",height = 25,width = 12)
#pheatmap(lmm_heatmap_plot,annotation_row = annotation_lmm,cluster_cols = F,clustering_distance_rows = "euclidean",
#annotation_colors = anno_col,cellheight = 6,fontsize_row = 5,fontsize_col = 10, border_color = NA)
#Heatmap(lmm_heatmap_plot_scale,clustering_distance_rows= "spearman",col = pheatmap_colors,cluster_columns = F,right_annotation = ha)#use pheatmap color
Heatmap(lmm_heatmap_plot_scale2,clustering_distance_rows= "spearman",cluster_columns = F,right_annotation = ha,name = "Model")
dev.off()

# quadratic heatmap
qua_heatmap_plot2<-eff_size[which(eff_size$sqtl %in% qua_sqtl2$sqtl),]
rownames(qua_heatmap_plot2)<-NULL
qua_heatmap_plot2<-column_to_rownames(qua_heatmap_plot2,var = "sqtl")
qua_heatmap_plot_scale2<-t(scale(t(qua_heatmap_plot2),scale = F))



qua_sqtl2<-qua_sqtl
qua_sqtl2$model[which(qua_sqtl2$sqtl %in% intersect)]<-"Linear & Quadratic"
annotation1<-rbind(lmm_sqtl2[,c("sqtl","model")],qua_sqtl2[,c("sqtl","model")])
rownames(annotation1)<-NULL
annotation1<-column_to_rownames(annotation1,var = "sqtl")


#pheatmap_colors <- rev(brewer.pal(n = 7, name = "RdYlBu"))
#pheatmap_colors <- colorRampPalette(pheatmap_colors)(300)





#plot1-qua
anno_col2<-list(model1=c("Linear & Quadratic"="#66a61e","quadratic model"="#8da0cb"),
                model2=c("quadratic model"="#8da0cb"))

#annotation for heatmap
c("Linear & Quadratic"="#8da0cb","linear model" = "#e78ac3","linear model"="#a6d854")
annotation_qua1<-qua_sqtl2[,c("sqtl","model")]
annotation_qua2<-qua_sqtl[,c("sqtl","model")]
annotation_qua<-merge(annotation_qua1,annotation_qua2,by="sqtl")
names(annotation_qua)<-c("sqtl","model1","model2")
annotation_qua<-column_to_rownames(annotation_qua,var="sqtl")
qua_heatmap_plot<-column_to_rownames(qua_heatmap_plot,var="sqtl")

#annotation for complexheatmap
col_anno<-c("Linear & Quadratic"="#66a61e","quadratic model"="#e7298a")
annotation_qua$color<-col_anno[annotation_qua$model1]
annotation_qua$sqtl <- rownames(annotation_qua)
ha = rowAnnotation(bar = annotation_qua$sqtl,
                   col = list(bar=setNames(as.character(annotation_qua$color), as.character(annotation_qua$sqtl))),
                   show_legend=c(bar=FALSE))



#qua_heatmap_plot_scale<-t(scale(t(qua_heatmap_plot),scale = F))

pdf("Fig4.qua.complexheatmap_scale_2.pdf",height = 6,width = 12)
#pheatmap(qua_heatmap_plot,annotation_row = annotation_qua,cluster_cols = F,clustering_distance_rows = "euclidean",
#         annotation_colors = anno_col2,cellheight = 6,fontsize_row = 5,fontsize_col = 10, border_color = NA)
Heatmap(qua_heatmap_plot_scale2,clustering_distance_rows= "spearman",col = pheatmap_colors,cluster_columns = F,right_annotation = ha)
Heatmap(qua_heatmap_plot_scale2,clustering_distance_rows= "spearman",cluster_columns = F,right_annotation = ha)
dev.off()

###############################################################################################################

## plot effect size box plot
# To use this function, you need to prepare some input files:
# phenotype data (do read as "phenotype")
# genotype data (do read as "genotype")
# and your sQTL file with p-value
#' x-genotype

library(ggpubr)
boxplotalleleeffect <- function(pheno,geno,model) {
  library(ggpubr)
  #prepare file
  geno1<-genotype[geno,]
  pheno1<-as.numeric(phenotype[pheno,])
  model=model
  file<-rbind(geno1,pheno1) %>% t() %>% as.data.frame()
  file$pseudotime<-str_split_fixed(rownames(file),"\\.",2)[,2]
  #change 0|1|2 into base type
  ref<-str_split_fixed(rownames(geno1),":",4)[,3]
  alt<-str_split_fixed(rownames(geno1),":",4)[,4]
  file$genotype<-"a"
  file$genotype[grep("0",file[,1])]<-paste0(ref,ref)
  file$genotype[grep("1",file[,1])]<-paste0(ref,alt)
  file$genotype[grep("2",file[,1])]<-paste0(alt,alt)
  file<-rownames_to_column(file,var = "sample")
  names(file)<-c("sample","geno","pheno","pseudotime","genotype")
  file$pheno<-as.numeric(file$pheno)
  file$genotype<-factor(file$genotype,levels = c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
  p1<-ggboxplot(file, x = "genotype", y ="pheno",
                color = "pseudotime", palette = "npg",font.x=c(9),font.y=c(9),xlab = geno,ylab=pheno) + 
    font("legend.title",size=8)+font("legend.text",size=8)+
    labs(title = model)
  return(p1)
}

library(ggplotify)
library(patchwork)

#only linear
p<-list()
for (i in 1:nrow(only_lmm)) {
  geno<-only_lmm$variants[i]
  pheno<-only_lmm$phenotye[i]
  p1<-boxplotalleleeffect(pheno = pheno,geno = geno,model="LMM-linear")
  p[[i]]<-as.ggplot(p1)
}
pdf("/Users/tongyihan/project/AIDA/data/figure/figure4/Fig4.LMM.linear.boxplot.pdf",width = 20,height = 96)
wrap_plots(p,ncol=4)
dev.off()

#linear & quadratic
p<-list()
for (i in 1:nrow(intersect_qtl)) {
  geno<-intersect_qtl$variants[i]
  pheno<-intersect_qtl$phenotye[i]
  p1<-boxplotalleleeffect(pheno = pheno,geno = geno,model="Linear & Quadratic")
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.Linear_Quadratic.boxplot.pdf",width = 20,height = 4)
wrap_plots(p,ncol=4)
dev.off()

#only quadratic
p<-list()
for (i in c(1:nrow(only_qua))) {
  geno<-only_qua$variants[i]
  pheno<-only_qua$phenotye[i]
  p1<-boxplotalleleeffect(pheno = pheno,geno = geno,model="Quadratic")
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.Quadratic.boxplot.pdf",width = 16,height = 8)
wrap_plots(p,ncol=4)
dev.off()

###############################################################################################################
#' boxplot by quantile
boxplotalleleeffect2 <- function(pheno,geno,model) {
  library(ggpubr)
  #prepare file
  geno1<-genotype[geno,]
  pheno1<-as.numeric(phenotype[pheno,])
  model=model
  file<-rbind(geno1,pheno1) %>% t() %>% as.data.frame()
  file$pseudotime<-str_split_fixed(rownames(file),"\\.",2)[,2]
  #change 0|1|2 into base type
  ref<-str_split_fixed(rownames(geno1),":",4)[,3]
  alt<-str_split_fixed(rownames(geno1),":",4)[,4]
  file$genotype<-"a"
  file$genotype[grep("0",file[,1])]<-paste0(ref,ref)
  file$genotype[grep("1",file[,1])]<-paste0(ref,alt)
  file$genotype[grep("2",file[,1])]<-paste0(alt,alt)
  file<-rownames_to_column(file,var = "sample")
  names(file)<-c("sample","geno","pheno","pseudotime","genotype")
  file$pheno<-as.numeric(file$pheno)
  file$genotype<-factor(file$genotype,levels = c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
  file$pseudotime<-as.numeric(gsub("Q","",file$pseudotime))
  colors<-c("#fb8072","#80b1d3","#fdb462")
  p1<-ggboxplot(file, x = "pseudotime", y ="pheno",
                fill = "genotype", palette = c("#fb8072","#80b1d3","#fdb462"),font.x=c(9),font.y=c(9),xlab = geno,ylab=pheno) + 
    font("legend.title",size=8)+font("legend.text",size=8)+
    labs(title = model)
  return(p1)
}

#only linear
p<-list()
for (i in 1:nrow(only_lmm)) {
  geno<-only_lmm$variants[i]
  pheno<-only_lmm$phenotye[i]
  p1<-boxplotalleleeffect2(pheno = pheno,geno = geno,model="LMM-linear")
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.linear.boxplot_v2.pdf",width = 20,height = 96)
wrap_plots(p,ncol=4)
dev.off()

#linear & quadratic
p<-list()
for (i in 1:nrow(intersect_qtl)) {
  geno<-intersect_qtl$variants[i]
  pheno<-intersect_qtl$phenotye[i]
  p1<-boxplotalleleeffect2(pheno = pheno,geno = geno,model="Linear & Quadratic")
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.Linear_Quadratic.boxplot_v2.pdf",width = 20,height = 4)
wrap_plots(p,ncol=4)
dev.off()

#only quadratic
p<-list()
for (i in c(1:nrow(only_qua))) {
  geno<-only_qua$variants[i]
  pheno<-only_qua$phenotye[i]
  p1<-boxplotalleleeffect2(pheno = pheno,geno = geno,model="Quadratic")
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.Quadratic.boxplot_v2.pdf",width = 16,height = 8)
wrap_plots(p,ncol=4)
dev.off()


#' dot plot 
#' input: raw perind count (after filtered sample you need and genotype that is significant)
library(tidyverse)
library(gganimate)
library(ggplot2)
library(reshape2)
#genotype
genotype<-read.table("merge_all_chr_unique_geno.tsv",header = T)
genotype<-genotype[,-c(1:2)]
genotype<-column_to_rownames(genotype,var = "ID")

#raw phenotype
#phenotype<-read.csv("/Users/tongyihan/project/AIDA/data/dynamic sQTL/test JP/JP_test_perind.counts.gz.qqnorm_chr22",sep="\t")
perind<-read.table("07.joint_all_perind.counts.only_2618.tsv",sep="\t",header=T)
#perind<-column_to_rownames(perind,var="clu")


genodotplot_df <- function(geno,pheno) {
  #prepare file
  ref<-str_split_fixed(geno,":",4)[,3]
  alt<-str_split_fixed(geno,":",4)[,4]
  geno1<-genotype[geno,]
  ref_ref<-names(geno1[,which(apply(geno1, 2, function(x) all(x == 0)))])
  ref_alt<-names(geno1[,which(apply(geno1, 2, function(x) all(x == 1)))])
  alt_alt<-names(geno1[,which(apply(geno1, 2, function(x) all(x == 2)))])
  perind1<-perind[which(perind$clu==pheno),]
  perind1<-melt(perind1,id.vars = "clu")
  perind1<-perind1[-which(perind1$value=="0/0"),] # delete sample with "0/0"
  perind1$numerator<-str_split_fixed(perind1$value,"/",2)[,1]
  perind1$dominator<-str_split_fixed(perind1$value,"/",2)[,2]
  perind1$quantile<-str_split_fixed(perind1$variable,"\\.",2)[,2]
  # extract perind value from certain type (ref_ref, ref_alt..)
  ref_ref_perind<-perind1[which(perind1$variable %in% ref_ref),]
  ref_ref_perind <-ref_ref_perind %>% group_by(quantile) %>% summarise(sum_numerator=sum(as.numeric(numerator)),sum_dominator=sum(as.numeric(dominator)))
  ref_ref_perind <- ref_ref_perind %>% mutate(names=paste0(ref,ref,"_",quantile),raw_count_ratio=sum_numerator/sum_dominator,sqtl=paste0(pheno,"_",geno))
  
  ref_alt_perind<-perind1[which(perind1$variable %in% ref_alt),]
  ref_alt_perind <-ref_alt_perind %>% group_by(quantile) %>% summarise(sum_numerator=sum(as.numeric(numerator)),sum_dominator=sum(as.numeric(dominator)))
  ref_alt_perind <- ref_alt_perind %>% mutate(names=paste0(ref,alt,"_",quantile),raw_count_ratio=sum_numerator/sum_dominator,sqtl=paste0(pheno,"_",geno))
  
  alt_alt_perind<-perind1[which(perind1$variable %in% alt_alt),]
  alt_alt_perind <-alt_alt_perind %>% group_by(quantile) %>% summarise(sum_numerator=sum(as.numeric(numerator)),sum_dominator=sum(as.numeric(dominator)))
  alt_alt_perind <- alt_alt_perind %>% mutate(names=paste0(alt,alt,"_",quantile),raw_count_ratio=sum_numerator/sum_dominator,sqtl=paste0(pheno,"_",geno))
  merge<-rbind(ref_ref_perind,ref_alt_perind,alt_alt_perind)
  return(merge)
}

color<-c("Q1"="#8dd3c7",
         "Q2" = "#ffffb3",
         "Q3" = "#bebada",
         "Q4" = "#fb8072",
         "Q5" = "#80b1d3",
         "Q6" = "#fdb462")
geno="chr21:34019120:G:C"
pheno="chr21:34073745:34074368_+"

test<-genodotplot_df(geno,pheno)
test$quantile<-factor(test$quantile,levels = c(paste0("Q",1:6)))
test$names<-factor(test$names,levels = rev(c(paste0(ref,ref,paste0("_Q",1:6)),paste0(ref,alt,paste0("_Q",1:6)),paste0(alt,alt,paste0("_Q",1:6)))))
ggplot(test,aes(x=sqtl,y=names,color=quantile,size=raw_count_ratio))+
  geom_point()+
  scale_color_manual(values = color)+
  labs(title = pheno)+
  theme_bw()

#only linear
p<-list()
for (i in 1:nrow(only_lmm)) {
  geno<-only_lmm$variants[i]
  pheno<-only_lmm$phenotye[i]
  pheno<-gsub("−","-",pheno)
  ref<-str_split_fixed(geno,":",4)[,3]
  alt<-str_split_fixed(geno,":",4)[,4]
  test<-genodotplot_df(geno,pheno)
  test$names<-factor(test$names,levels = rev(c(paste0(ref,ref,paste0("_Q",1:6)),paste0(ref,alt,paste0("_Q",1:6)),paste0(alt,alt,paste0("_Q",1:6)))))
  p1<-ggplot(test,aes(x=sqtl,y=names,color=quantile,size=raw_count_ratio))+
    geom_point()+
    scale_color_manual(values = color)+
    labs(title = pheno)+
    theme_bw()
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.linear.dotplot2.pdf",width = 20,height = 96)
wrap_plots(p,ncol=4)
dev.off()

#linear & quadratic
p<-list()
for (i in 1:nrow(intersect_qtl)) {
  geno<-intersect_qtl$variants[i]
  pheno<-intersect_qtl$phenotye[i]
  pheno<-gsub("−","-",pheno)
  ref<-str_split_fixed(geno,":",4)[,3]
  alt<-str_split_fixed(geno,":",4)[,4]
  test<-genodotplot_df(geno,pheno)
  test$names<-factor(test$names,levels = rev(c(paste0(ref,ref,paste0("_Q",1:6)),paste0(ref,alt,paste0("_Q",1:6)),paste0(alt,alt,paste0("_Q",1:6)))))
  p1<-ggplot(test,aes(x=sqtl,y=names,color=quantile,size=raw_count_ratio))+
    geom_point()+
    scale_color_manual(values = color)+
    labs(title = pheno)+
    theme_bw()
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.Linear_Quadratic.dotplot.pdf",width = 20,height = 4)
wrap_plots(p,ncol=4)
dev.off()

#only quadratic
p<-list()
for (i in c(1:nrow(only_qua))) {
  geno<-only_qua$variants[i]
  pheno<-only_qua$phenotye[i]
  ref<-str_split_fixed(geno,":",4)[,3]
  alt<-str_split_fixed(geno,":",4)[,4]
  test<-genodotplot_df(geno,pheno)
  test$names<-factor(test$names,levels = rev(c(paste0(ref,ref,paste0("_Q",1:6)),paste0(ref,alt,paste0("_Q",1:6)),paste0(alt,alt,paste0("_Q",1:6)))))
  p1<-ggplot(test,aes(x=sqtl,y=names,color=quantile,size=raw_count_ratio))+
    geom_point()+
    scale_color_manual(values = color)+
    labs(title = pheno)+
    theme_bw()
  p[[i]]<-as.ggplot(p1)
}
pdf("Fig4.LMM.Quadratic.dotplot.pdf",width = 16,height = 8)
wrap_plots(p,ncol=4)
dev.off()










