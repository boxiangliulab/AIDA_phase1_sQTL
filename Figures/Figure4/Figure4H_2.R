#   ____________________________________________________________________________
#   Script information                                                      ####

# title: plot significant chnage of introns usages
# author: Yihan Tong
# date: 2023-03-09
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# conda activate r-4.2.3

#' this script is for plot significant chnage of introns usages
#' input: Bcell_quantile_only_pc_and_lncRNA_perind.counts + all_sample_ANOVA_value.csv

library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(reshape2)
library(pheatmap)
library(ggplotify)
library(patchwork)


perind_test<-read.csv("introns_from_perind_counts_in_qqnorm.tsv",sep="\t")
anova_df<-read.csv("all_sample_ANOVA_value.csv")
anova_df<-anova_df[,-1]
anova_df$p_adj<-p.adjust(anova_df$p_value)

# plot p.adj distribution
anova_df$rank<-rank(anova_df$p_adj)

ggplot(anova_df,aes(x=rank,y=p_adj))+geom_point()+geom_hline(yintercept = 1,color="red",linetype="dashed")+theme_classic()
ggsave("01.all_sample_ANOVA_pvalue_point.pdf",width = 4,height = 4)

ggplot(anova_df,aes(x=p_adj))+geom_density()+theme_classic()+geom_vline(xintercept = 0.05,linetype="dashed",color="red")
ggsave("01.all_sample_ANOVA_pvalue_distribution.pdf",width = 4,height = 4)

print("finsh p-adj value plot")

# filter p.adj <= 0.05 #

anova_df<-dplyr::filter(anova_df,p_adj<=0.05)
print(paste0("There are ", nrow(anova_df)," introns whose p.adj is smaller than 0.05"))


index<-which(perind_test$chrom %in% anova_df$junction)
perind_test_extr<-perind_test[index,]
rownames(perind_test_extr)<-NULL
perind_test_extr2<-melt(perind_test_extr,id.vars="chrom")

perind_test_extr2<-dplyr::filter(perind_test_extr2,value!="0/0")
#conver fraction to decimal#
perind_test_extr2$ratio <- sapply(strsplit(perind_test_extr2$value, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

#add quantile info#
perind_test_extr2$quantile<-str_split_fixed(perind_test_extr2$variable,"\\.",5)[,2]

#calculate mean of each intron junction in all the quantiles it has
data<-data.frame()
for (i in unique(perind_test_extr2$chrom)) {
  df_tmp<-dplyr::filter(perind_test_extr2,chrom==i)
  for (q in unique(perind_test_extr2$quantile)) {
    data_tmp<-data.frame(intron=i,quantile=q,mean=mean(df_tmp$ratio[which(df_tmp$quantile==q)]),p_adj=anova_df$p_adj[which(anova_df$junction==i)])
    data<-rbind(data,data_tmp)
  }
}

write.table(data,"01.sig_intron_mean_usage_across_quantile.tsv",sep="\t",row.names = F,quote = F)
data2<-dcast(data=data,formula=intron~quantile,value.var="mean") ##column: quantile  row:introns
data2<-column_to_rownames(data2,var="intron")
write.table(data2,"01.sig_intron_mean_usage_across_quantile_heatmap.tsv",sep="\t",row.names = T,quote = F)
dev.off()

library(ComplexHeatmap)
scaled_mat = t(scale(t(data2)))
pdf("02.sig_intron_heatmap.pdf",height = 10,width = 6)
Heatmap(scaled_mat,clustering_distance_rows= "spearman")
dev.off()

#dev.off()

# plot boxplot

unique_intron<-unique(perind_test_extr2$chrom)
p<-list()
colors<-rev(c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
for (i in unique_intron) {
  df_tmp<-dplyr::filter(perind_test_extr2,chrom==i)
  num<-length(unique(df_tmp$quantile))
  colors2<-colors[1:num]
  m<-ggplot(df_tmp, aes(x = quantile, y = ratio)) +
    geom_boxplot(fill = colors2, color = "black", outlier.color = "black") +
    geom_point(size = 2, alpha = 0.2) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(size = 0.5, linetype = 'dashed'),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.line = element_line(size = 1, color = "black"),
      legend.position = "none",
      plot.title = element_text(size = 14, color = "black", face = "bold"),
      plot.subtitle = element_text(size = 12, color = "black"),
      plot.caption = element_text(size = 10, color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      strip.background = element_rect(fill = "gray90", color = NA),
      strip.text = element_text(size = 10, color = "black", face = "bold")
    ) +
    labs(
      x = "Quantile",
      y = "Intron Usage",
      title = "Boxplot of Intron Usage by Quantile",
      subtitle = i
    )
  median<-data.frame()
  for (q in unique(df_tmp$quantile)) {
    m_tmp<-dplyr::filter(df_tmp,quantile==q)
    median_tmp<-data.frame(median=median(m_tmp$ratio),quantile=q)
    median<-rbind(median,median_tmp)
  }
  median$quantile2<-as.numeric(str_split_fixed(median$quantile,"Q",2)[,2])
  median2<-median
  p1<-m + geom_point(data = median2, aes(x = quantile2, y = median), size = 4, color = "red")+
    geom_smooth(data = median2, aes(x = quantile2, y = median), method = "loess", color = "grey0", size = 2, se = FALSE)
  p[[i]]<-as.ggplot(p1)
}

pdf("03.intron_usage_boxplot.pdf",width = 50,height = 400)
wrap_plots(p,ncol=10)
dev.off()

# plot dot plot
target_introns<-read.table("introns_needed_to_show.txt")
perind_test_extr3<-melt(perind_test,id.vars="chrom")
perind_test_extr3<-dplyr::filter(perind_test_extr3,value!="0/0")
#分数转化成小数#
perind_test_extr3$ratio <- sapply(strsplit(perind_test_extr3$value, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

#add quantile info#
perind_test_extr3$quantile<-str_split_fixed(perind_test_extr3$variable,"\\.",5)[,2]
perind_test_extr3$clu<-str_split_fixed(perind_test_extr3$chrom,":",4)[,4]

p<-list()
colors<-rev(c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
p<-list()
colors<-rev(c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"))
for (i in target_introns$V1) {
  df<-dplyr::filter(perind_test_extr3,clu==i)
  for (j in unique(df$chrom)) {
    df_tmp<-dplyr::filter(df,chrom==j)
    num<-length(unique(df_tmp$quantile))
    colors2<-colors[1:num]
    mean<-data.frame()
    for (q in unique(df_tmp$quantile)) {
      m_tmp<-dplyr::filter(df_tmp,quantile==q)
      mean_tmp<-data.frame(mean=mean(m_tmp$ratio),quantile=q)
      mean<-rbind(mean,mean_tmp)
    }
    mean$intron=j
    mean$quantile<-factor(mean$quantile,levels = rev(c(paste0("Q",1:6))))
    #median$quantile2<-as.numeric(str_split_fixed(median$quantile,"Q",2)[,2])
    #median2<-median
    p1<-ggplot(mean,aes(x=intron,y=quantile,size=mean,color=quantile))+
      geom_point()+
      scale_color_manual(values = colors)+
      theme_classic() +
      labs(
        x = "Intron",
        y = "Quantile",
        title = i,
        subtitle = j
      )
    p[[j]]<-as.ggplot(p1)
  }
}

pdf("Fig4.intron_usage.dot_plot.pdf",width = 16,height = 160)
wrap_plots(p,ncol=4)
dev.off()




