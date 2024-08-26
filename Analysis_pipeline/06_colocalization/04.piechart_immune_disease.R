###This script was used to generate Figure7C & Figure7D. Figure7C shows examples of sGenes performing colocalization with five immune diseases Asthma, SLE, RA, AD, GD.
###Figure7D shows TCHP's gene expression, junction usage, sQTL, eQTL and H4PP conditions among 19 cell types.

###Figure7C
##plot star plot

library(graphics)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
data(diamonds)
dstar<- sample_n(diamonds, 100)
dstar$log10carat <- log10(dstar$carat)
dstar$log10price <- log10(dstar$price)
dstar<-dstar[order(dstar$cut,decreasing=T),]

stars(dstar[,2:6], key.loc = c(-2, 10), scale = TRUE, 
      locations = NULL, len =1, radius = TRUE,
      full = TRUE, labels = NULL,draw.segments = TRUE,
      col.segments=palette(brewer.pal(7,"Set1"))[1:5])

loc <- data.matrix(dstar[,11:12])
stars(dstar[,2:6], key.loc = c(-1, 3), scale = TRUE, 
      locations = loc, len =0.07, radius = TRUE,
      full = TRUE, labels = NULL, draw.segments = TRUE,
      col.segments=palette(brewer.pal(7,"Set1"))[1:5],
      frame.plot=TRUE,axes = TRUE, 
      xlab="log10(carat)", ylab="log10price",
      xlim=c(-0.7,0.7))


setwd("/Users/tongyihan/project/AIDA/data/figure/figure7/figure7B/")
clean_data <- function(file,disease) {
  file$gene<-rownames(file)
  file<-melt(file,id.vars = "gene") %>%  dplyr::filter(!is.na(value))
  names(file)<-c("gene","celltype","H4")
  file$disease<-disease
  return(file)
}

# SLE
sle<-read.table("SLE.txt",check.names = F,row.names = 1)
#sle$gene<-rownames(sle)
#sle<-melt(sle,id.vars = "gene") %>%  dplyr::filter(!is.na(value))
#names(sle)<-c("gene","celltype","H4")
#sle$disease<-"SLE"
sle<-clean_data(sle,disease = "SLE")
#RA
ra<-read.table("RA.txt",check.names = F,row.names = 1)
ra<-clean_data(ra,disease = "RA")
#GD
gd<-read.table("GD.txt",check.names = F,row.names = 1)
gd<-clean_data(gd,disease = "GD")
#AD
ad<-read.table("AD.txt",check.names = F,row.names = 1)
ad<-clean_data(ad,disease = "AD")
#Asthma
asthma<-read.table("Asthma.txt",check.names = F,row.names = 1)
asthma<-clean_data(asthma,disease = "Asthma")

merge<-rbind(sle,ra,gd,ad,asthma)
new_merge<-separate_rows(merge,gene,sep = ",")
new_merge <- new_merge[!(new_merge$gene %in% c("AC007834.1", "AL159163.1")), ]

# add TCHP to GD in naive CD4+ T cell
new_merge[75,1] <- "TCHP"
new_merge[75,2]<-"CD4+_T_naive"
new_merge[75,3]<-1
new_merge[75,4]<-"GD" 
# add TOP3B to GD in naive CD4+ T cell
new_merge[76,1] <- "TOP3B"
new_merge[76,2]<-"cDC2"
new_merge[76,3]<-1
new_merge[76,4]<-"SLE" 





#set coordinate

gene_order<-c(unique(new_merge$gene))
disease_order<-c("SLE","RA","GD","AD","Asthma")




for (i in seq_along(label_dict)) {
  new_merge$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],new_merge$celltype)
}

gene_vec <- 1:length(unique(new_merge$gene))
celltype_vec <- 1:length(unique(new_merge$celltype))






# here you can change the order of genes/celltypes
names(gene_vec) <- unique(new_merge$gene)
names(celltype_vec) <- unique(new_merge$celltype)

new_merge$gene_coord <- gene_vec[new_merge$gene]
new_merge$celltype_coord <- celltype_vec[as.character(new_merge$celltype)]
new_merge$pair<-paste0(new_merge$gene,"-",new_merge$celltype)
pair<-unique(new_merge$pair)
df<-data.frame(new_merge[,c(5:7,1,2)]) %>% unique()
for (i in disease_order) {
  df_tmp<-dplyr::filter(new_merge,disease==i)
  names(df_tmp)<-c("gene","celltype",i,"disease","gene_coord","celltype_coord","pair")
  diff_pair<-setdiff(pair,unique(df_tmp$pair))
  diff_pair_df<-data.frame()
  for (j in 1:length(diff_pair)) {
    diff_pair_tmp<-data.frame(gene=str_split_fixed(diff_pair[j],"-",2)[,1],
                              celltype=str_split_fixed(diff_pair[j],"-",2)[,2],
                              H4=0,
                              disease=i,
                              gene_coord=new_merge$gene_coord[which(new_merge$pair==diff_pair[j])[1]],
                              celltype_coord=new_merge$celltype_coord[which(new_merge$pair==diff_pair[j])[1]],
                              pair=diff_pair[j])
    diff_pair_df<-rbind(diff_pair_df,diff_pair_tmp) 
  }
  names(diff_pair_df)<-c("gene","celltype",i,"disease","gene_coord","celltype_coord","pair")
  df_tmp<-rbind(diff_pair_df,df_tmp)
  print(paste0(i," has ", nrow(df_tmp)," pairs"))
  df_tmp<-data.frame(df_tmp[,c(3,7)])
  names(df_tmp)<-c(i,"pair")
  df<-merge(df,df_tmp,by="pair")
}

df[, 6:10][df[, 6:10] != 0] <- 1
for (i in 1:nrow(df)) {
  df$SLE_pro[i]<-df$SLE[i]/sum(df[i,6:10])
  df$RA_pro[i]<-df$RA[i]/sum(df[i,6:10])
  df$GD_pro[i]<-df$GD[i]/sum(df[i,6:10])
  df$AD_pro[i]<-df$AD[i]/sum(df[i,6:10])
  df$Asthma_pro[i]<-df$AD[i]/sum(df[i,6:10])
}


df_long<-new_merge[,c(4:6)]

#new_merge$gene<-factor(new_merge$gene,levels = gene_order)
#new_merge$disease<-factor(new_merge$disease,levels = disease_order)

#loc <- data.matrix(df[,2:3])
#stars(df[,c(6:10)], key.loc = c(8:3), scale = FALSE,
 #     locations = loc, len =0.7, radius = TRUE,
  #    full = TRUE, labels = NULL, draw.segments = TRUE,
   #   col.segments=disease_col,
    #  frame.plot=TRUE,axes = FALSE, 
     # xlab="log10(carat)", ylab="log10price",
      # xlim=c(0,30))
#axis(side = 1, at = 0:10)
#axis(side = 2, at = 0:20, labels = paste0("a", 0:20))


###complex piechart

library(scatterpie)
library(ggforce)
disease_col<-c("AD"="#fc8d62","Asthma"="#66c2a5","GD"="#8da0cb","RA"="#e78ac3","SLE"="#a6d854")

df_long.new <- dcast(df_long, formula = gene_coord + celltype_coord ~ disease, fun.aggregate = length)
#df_long.new[74,] <- c(27,19,0,0,0,0,1) # add TOB3B in SLE of cDC2
df_long.new <- df_long.new %>% mutate_at(tail(names(df_long.new),5),as.numeric)
#gene_vec<-c(gene_vec,"TOP3B"=27)
#celltype_vec<-c(celltype_vec,"cDC2"=19)

ggplot() +
  geom_scatterpie(aes(x = celltype_coord, y = gene_coord,alpha=1,size=100), data = df_long.new, cols = c("AD", "Asthma", "GD", "RA", "SLE"), pie_scale = 1.2,size=2,color=NA) +
  coord_fixed() +
  scale_x_continuous(name = "Cell type",breaks = seq(1,max(df_long.new$celltype_coord),1),labels = names(celltype_vec)) +
  scale_y_continuous(name = "Gene",breaks = seq(1,max(df_long.new$gene_coord),1),labels = names(gene_vec)) +
  theme_bw() +
  theme(
    axis.line = element_line(size = 0.5, color = "black"),  
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14, face = "bold"),  
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45,hjust = 1)
  )+
  scale_fill_manual(values = c("AD"="#fc8d62","Asthma"="#66c2a5","GD"="#8da0cb","RA"="#e78ac3","SLE"="#a6d854"))+
  geom_circle(aes(x0 = celltype_coord, y0 = gene_coord, r = 0.45), data = df_long.new,  fill = NA, color = "black",inherit.aes = FALSE) +
  coord_equal()





ggsave("figure7B.piechart.v2.pdf",width = 8,height = 10)



### bar chart
#new_merge[82,]<-c("TOP3B","cDC2","unkown","SLE",27,19,"TOP3B-cDC2")

new_merge2<-new_merge[,c(4,1)] %>% unique()
bar_plt<-as.data.frame(table(new_merge2$disease))
bar_plt$Var1<-factor(bar_plt$Var1,levels = c("AD","GD","Asthma","RA","SLE"))
Freq<-c(5,10,5,16,17)
Var1<-c("AD","Asthma","GD","RA","SLE")
bar_plt<-data.frame(Var1=Var1,Freq=Freq)
ggplot(bar_plt,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("AD"="#fc8d62","Asthma"="#66c2a5","GD"="#8da0cb","RA"="#e78ac3","SLE"="#a6d854"))+
  theme_classic()+guides(fill=FALSE)
ggsave("figure7B.barplot.v2.pdf",width = 6,height = 8)




write.table(df_long.new,"figure7B.plot.df.v2.tsv",sep="\t",row.names = F)
write.table(new_merge,"merge_all_disease_ori.v2.tsv",sep="\t",row.names = F)



###save csv to pdf
library(gridExtra)
sum_table<-new_merge %>% group_by(gene,disease) %>% summarise(celltype_num=length(celltype))
sum_table_wide<-dcast(sum_table,formula = gene ~ disease)
sum_table_wide[is.na(sum_table_wide)] <- 0
gene_order<-rev(c("ELF1","TPCN2","IKBKE","SIPA1","IRF5","ARHGAP15","THAP7","C11orf80","CD40","AFF3",
              "AP003774.2","RNASET2","AL159163.1","HDGF","ARID5B","PIGT","CD83","FCRL3",
              "CCT7","RER1","AC007834.1","TCHP","CMTM7","JMJD1C","STAT6","ORMDL3","TOP3B"))
sum_table_wide<-sum_table_wide[match(gene_order,sum_table_wide$gene),]


#sum_table_wide$gene<-factor(sum_table_wide$gene,levels = c("ELF1","TPCN2","IKBKE","SIPA1","IRF5","ARHGAP15","THAP7","C11orf80","CD40","AFF3",
 #                                                     "AP003774.2","RNASET2","AL159163.1","HDGF","ARID5B","PIGT","CD83","FCRL3",
  #                                                    "CCT7","RER1","AC007834.1","TCHP","CMTM7","JMJD1C","STAT6","ORMDL3","TOP3B"))


sum_table_wide_new <- tableGrob(sum_table_wide, rows = NULL)


out<-arrangeGrob(sum_table_wide_new,respect = TRUE)
out2 <- grid.arrange(out,heights=c(0.6,0.4))
print(out2)
ggsave('sum_table.pdf', width = 4, height = 18)


write.csv(sum_table_wide,"sum_table_wide.csv",row.names = F)


###Figure7D
##plot TCHP junction usage
cellnum<-read.table("C:/users/90410/desktop/aida/supple_file/SIPA1/cell_type_num_per_individual.txt",sep=" ",header=TRUE)
junc<-read.table("C:/users/90410/desktop/aida/TCHP/junction_usage_TCHP.txt",sep="\t")
stat_junc<-as.data.frame(matrix(NA,19,3))
rownames(stat_junc)<-celltype_random
cellnum$V4<-paste(cellnum$DCP_ID,cellnum$Annotation_Level2,"bam",sep=".")
rownames(cellnum)<-cellnum$V4
for(i in 1:19){
  num<-which(grepl(celltype_random[i],junc[,1],fixed=TRUE)==TRUE)
  pseudo<-length(num)
  total_read<-sum(junc[num,2])
  stat_junc[i,1]<-total_read
  stat_junc[i,2]<-pseudo
  stat_junc[i,3]<-total_read/pseudo
}
color = colorRampPalette(c("navy", "white", "red"))(50)
stat_junc[,1]<-stat_junc[,3]
stat_junc[,2]<-stat_junc[,3]
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC",
                           "#053061"))(300) %>% rev()
pheatmap(as.data.frame(stat_junc[,1]),cluster_cols = FALSE,cluster_rows = FALSE, color =col2,cellwidth=9,cellheight = 9,
         filename = 'test.pdf', width = 3, height = 8,show_rownames = F,show_colnames = F)
###plot TCHP gene expression
gene<-read.table("C:/users/90410/desktop/aida/TCHP/zhangyuntian_TCHP_19_celltype_mean_expression.tsv",sep="\t",header=TRUE)
pheatmap(as.data.frame(gene[,2]),cluster_cols = FALSE,cluster_rows = FALSE, color =col2,cellwidth=9,cellheight = 9,
         filename = 'test.pdf', width = 3, height = 8,show_rownames = F,show_colnames = F)
