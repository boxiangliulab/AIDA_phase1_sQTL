##manhattan plot CLEC2D
##compare original female & now female
female_origin<-read.table("/Users/zhangyuntian/Desktop/aida_remake/aida/CLEC2D/male_female_total_12_12/original_CD4+_T_em_female_12_CLEC2D.txt",
                          sep=" ")
female_now<-read.table("/Users/zhangyuntian/Desktop/aida_remake/aida/CLEC2D/male_female_total_12_12/CD4+_T_em_female_CLEC2D.txt",
                 sep=" ")
plot(female_origin$V12,female_now$V12,xlab="combine_CLEC2D_chr12_9670678_9687902_pval",
     ylab="seperate_CLEC2D_chr12_9670678_9687902_pval")
plot(female_origin$V13,female_now$V13,xlab="combine_CLEC2D_chr12_9670678_9687902_beta",
     ylab="seperate_CLEC2D_chr12_9670678_9687902_beta")
library(ggplot2)
library(gcookbook)
library(dplyr)
GWAS<-read.table("/Users/zhangyuntian/Desktop/aida_remake/aida/CLEC2D/male_female_total_12_12/CD4+_T_em_female_CLEC2D.txt",sep=" ")
GWAS<-read.table("/Users/zhangyuntian/Desktop/aida_remake/aida/CLEC2D/GWAS_locuszoom")
newGWAS<-as.data.frame(matrix(NA,nrow(GWAS),5))
for(i in 1:nrow(GWAS)){
  newGWAS[i,1]<-strsplit(GWAS[i,2],"r")[[1]][2]
  newGWAS[i,2]<-strsplit(GWAS[i,8],":")[[1]][2]
  newGWAS[i,3]<-strsplit(GWAS[i,8],":")[[1]][3]
  newGWAS[i,4]<-strsplit(GWAS[i,8],":")[[1]][4]
  newGWAS[i,5]<-GWAS[i,12]
}
GWAS<-newGWAS
rsid<-read.table("/Users/zhangyuntian/Desktop/aida_remake/aida/CLEC2D/CLEC2D_snp_anno.txt",sep="\t",header=TRUE)
GWAS$V1<-paste(GWAS$V1,GWAS$V2,sep=":")
rsid$V1<-paste(rsid$Chromosome,rsid$Position,sep=":")


total<-merge(GWAS,rsid,by="V1")
total$LDscore<-0
##total<-total[,-2]
LD<-read.table("/Users/zhangyuntian/Desktop/aida_remake/aida/CLEC2D/LDmatrix.txt",sep="\t",header=TRUE)
num<-which(LD[,1]=="rs3764022")

LD<-LD[num,]
for(i in 2:ncol(LD)){
  num<-which(total[,7]==colnames(LD)[i])
  if(length(num)){
    total[num,17]<-LD[1,i]}
}


color<-c("#EE4000","#EEAD0E","#33cc00","#00ccff","#6959CD")
total<-na.omit(total)
for(i in 1:nrow(total)){
  if((total[i,17]>=0)&&(total[i,17]<0.2))total[i,1]<-color[5]
  if((total[i,17]>=0.2)&&(total[i,17]<0.4))total[i,1]<-color[4]
  if((total[i,17]>=0.4)&&(total[i,17]<0.6))total[i,1]<-color[3]
  if((total[i,17]>=0.6)&&(total[i,17]<0.8))total[i,1]<-color[2]
  if((total[i,17]>=0.8)&&(total[i,17]<=1))total[i,1]<-color[1]
}
total$V1<-factor(total$V1,levels=color)
total$Position<-total$Position/1000000
num<-which(total[,7]=="rs3764022")
p1<-ggplot(total,aes(x=Position,y=-log(V5,10),fill=V1))+geom_point(size=2.5,shape=21)+theme_classic()+
  scale_fill_manual(values=c("#EE4000","#EEAD0E","#33cc00","#00ccff","#6959CD"))+xlim(9.6,10)+ylim(0,25)+
  xlab(NULL)+ylab("-log10(P)")+theme(panel.border = element_rect(color = "black", size = 1, fill = NA))+
  ggrepel::geom_text_repel(aes(label=dbSNP,color=V1),total[num,],
                           size = 4, #ע???ı?????????С
                           box.padding = 0.5, #?ֵ????ľ???
                           point.padding = 0, #?ֵ????ľ??룬????Χ?Ŀհ׿???
                           min.segment.length = 0, #???߶ο???ʡ??
                           segment.color = "black", #segment.colour = NA, ????ʾ?߶?
                           show.legend = F)+guides(fill=FALSE)

p<-ggarrange(p1,p4,p2,p3,nrow=4,ncol=1)
ggsave("coloc_plot_new.pdf",p,width=5,height=10)



##read in correlation between GWAS sample size and colocalization proportion

library(tidyverse)
library(data.table)
library(fs)
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggpmisc)
library(ggplot2)
dat<-read_xlsx("C:/users/90410/desktop/aida/supple_file/corr_sample_proportion_coloc.xlsx")
formula_2 <- y ~ x
p2 <- ggplot(dat, aes(x = size, y = proportion)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method="lm",formula=formula_2)+
  stat_poly_eq(formula = formula_2,parse=TRUE) +
  xlab("sample size") +
  ylab("proportion")

ggsave("C:/users/90410/desktop/aida/supple_file/coloc_size.pdf", p2, height = 3, width = 3)
