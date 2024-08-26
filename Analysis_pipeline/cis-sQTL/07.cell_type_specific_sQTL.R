###This script was used to plot Figure3J and select cell-type specific sQTLs related to SLE (lfsr<0.05).
library(readxl)
snp<-read_xlsx("C:/users/90410/desktop/aida/supple_file/data/gene_snp_lfsr.xlsx")
gene<-read_xlsx("C:/users/90410/desktop/aida/supple_file/table/SLE_related_gene.xlsx")
common<-intersect(gene$gene,snp$gene)
SLE<-c("IL7R","IL16","IL32",                ##interleukin
       "IRF5","IRF7","TMEM173","MX1",       ##interferon-related
       "CCL4","CCL5",                       ##chemokine
       "ATG5","BTG3","PTTG1","CARD8",       ##cell cycle related
       "ARL5B","STAT4","LTA","FAS","CD40","CD83",    ##signalling transduction
       "ELF1","XBP1","BACH2","ARID5B",      ##transcription factor
       "AK3","ADARB1","CAT","TPMT","UBASH3A","TRIM21","TIMMDC1")      ##other genes related to immunology
annotation_row = data.frame(
  GeneClass = factor(rep(c("Immune related","Cell cycle related","Signal transduction",
                           "Transcription factors","Other genes"), c(9,4, 6,4,7)))
)
annotation_colour=c("#fbb4ae","#decbe4","#fed9a6","#ffffcc","#e5d8bd")
names(annotation_colour)<-c("Immune related","Cell cycle related","Signal transduction",
                            "Transcription factors","Other genes")
celltype_random=c("CD14+ Monocyte","CD16+ Monocyte","CD16+ NK","cm CD4+ T","cyt CD4+ T","em CD4+ T","Naive CD4+ T",
                  "CD56+ NK","GZMB+ CD8+ T","GZMK+ CD8+ T","Naive CD8+ T","IGHMhi memory B","IGHMlo memory B","MAIT","Treg",
                  "cDC2","GZMK+ gdT","GZMK- gdT","Naive B")
rownames(snp)<-snp$gene
snp<-as.data.frame(snp)
select.snp<-as.data.frame(snp[SLE,])

select.snp<-select.snp[,-c(1,2)]
rownames(select.snp)<-SLE
colnames(select.snp)<-celltype_random
for(i in 1:length(SLE)){
  for(j in 1:19){
    new<-select.snp[i,j]
    if(new<0.05)select.snp[i,j]<-1
    if(new>=0.05)select.snp[i,j]<-0
  }
}
library(pheatmap)
library(ggplot2)
colours = colorRampPalette(c("#F5F5F5", "#4682B4"))(2)
rownames(annotation_row) = rownames(select.snp)
ann_colors = list(
  GeneClass = c(`Immune related` = "#fbb4ae",`Cell cycle related`="#decbe4",
                `Transcription factors`="#ffffcc",
                `Signal transduction`="#fed9a6",`Other genes`="#d3d3d3")
)
p<-pheatmap(select.snp,border_color = "grey",color=colours,cluster_row=FALSE,angle_col=45,annotation_row = annotation_row,
            annotation_colors = ann_colors,
         cluster_col=FALSE,cellwidth = 9,cellheight = 9,legend=FALSE,fontsize_row=10,fontsize_col = 8)
ggsave("C:/users/90410/desktop/aida/supple_file/picture/figure3/SLE_related_gene_update_7.14.pdf",p,width=8,height=6)

###plot the violin plot of cell-type specific sQTLs in CD83.

intron_new=c("chr6%3A14118065%3A14124955%3Aclu_30543_+")
SNP_new=c("chr6_14127323.txt")
i<-1
j<-1
for(i in 1:length(intron)){
  for(j in 1:length(celltype)){
    intron<-read.table(paste("C:/Users/90410/Desktop/aida/supple_file/CD83/",celltype[j],intron_new[i],"_phe.txt",sep=""),sep="\t")
    col_intron<-read.table(paste("C:/Users/90410/Desktop/aida/sbQTL/phename/",celltype[j],"_phename.txt",sep=""),sep=" ")
    colnames(intron)<-col_intron[1,]
    SNP<-read.table(paste("C:/Users/90410/Desktop/aida/supple_file/CD83/",SNP_new[i],sep=""),sep="\t")
    col_SNP<-read.table("C:/Users/90410/Desktop/aida/sbQTL/phename/geno.txt",sep="\t")
    colnames(SNP)<-col_SNP[1,]
    inter<-intersect(colnames(intron),colnames(SNP))
    geno<-SNP[,inter]
    phe<-intron[,inter]
    total<-t(rbind(geno,phe))
    ref<-"T"
    alt<-"C"
    for(k in 1:nrow(total)){
      if(startsWith(total[k,1],'0|0'))total[k,1]<-paste(ref,ref,sep="")
      if(startsWith(total[k,1],'0|1'))total[k,1]<-paste(ref,alt,sep="")
      if(startsWith(total[k,1],'1|0'))total[k,1]<-paste(ref,alt,sep="")
      if(startsWith(total[k,1],'1|1'))total[k,1]<-paste(alt,alt,sep="")
    }
    total<-as.data.frame(total)
    colnames(total)<-c("geno","pheno")
    total$geno<-factor(total$geno,levels=c(paste(ref,ref,sep=""),paste(ref,alt,sep=""),paste(alt,alt,sep="")))
    total<-total[-1,]
    total$pheno<-as.numeric(total$pheno)
  }
}
color=c("#fee090","#fdae61", "#bf812d",
        "#9e9ac8","#6a51a3","#807dba",
        "#bcbddc","#8c510a","#3690c0","#74a9cf",
        "#a6bddb", "#a1d99b","#74c476",
        "#0570b0","#CCCC4D","#f46d43",
        "#f768a1", "#dd3497","#c7e9c0")
p1<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[12],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("IGHMhi memory B")   ##(p = 0.54)

p2<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[13],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("IGHMlo memory B") ##(p = 2.6e-06)
  
p3<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[19],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("Naive B")  ##(p = 0.86)

p4<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[7],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("Naive CD4+ T")  ## (p = 0.58)

p5<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[11],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("Naive CD8+ T")  ##(p = 0.81)

p6<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[16],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("cDC2")  ##(p = 0.5)

p7<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[17],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("GZMK+ gdT")  ##(p = 0.76)

p8<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[1],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("CD14+ Monocyte") ##(p = 0.22)

p9<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[2],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("CD16+ Monocyte")
p10<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[4],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("cm CD4+ T")
p11<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[6],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("em CD4+ T")
p12<-ggplot(total,aes(x=geno,y=pheno))+
  geom_violin(aes(color=geno),trim=FALSE)+
  geom_boxplot(aes(color=geno),width=0.05)+
  scale_color_manual(values=c(rep(color[3],3)))+theme_classic()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_text(size=14,hjust=0.5),
        axis.text = element_text(size=14))+
  ggtitle("CD16+ NK")

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=4,nrow=3)

ggsave("CD83_twelve.pdf",width=6,height=4.5)
