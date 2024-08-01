library(factoextra)
library(ggplot2)
library(reshape2)
library(corrplot)

##change the cell name to pretty names
#' code to conver cell type names

label_dict<-list("CD14\\+_Monocyte"="CD14+ Monocyte",
                 "CD16\\+_Monocyte"="CD16+ Monocyte",
                 "cDC2"="cDC2",
                 "CD16\\+_NK"="CD16+ NK",
                 "CD56\\+_NK"="CD56+ NK",
                 "CD4\\+_T_cm"="cm CD4+ T",
                 "CD4\\+_T_cyt"="cyt CD4+ T",
                 "CD4\\+_T_em"="em CD4+ T",
                 "CD4\\+_T_naive"="Naive CD4+ T",
                 "Treg"="Treg",
                 "CD8\\+_T_naive"="Naive CD8+ T",
                 "CD8\\+_T_GZMB\\+"="GZMBhi CD8+ T",
                 "CD8\\+_T_GZMK\\+"="GZMKhi CD8+ T",
                 "MAIT"="MAIT",
                 "gdT_GZMK\\-"="GZMBhi gdT",
                 "gdT_GZMK\\+"="GZMKhi gdT",
                 "naive_B"="Naive B",
                 "IGHMlo_memory_B"="IGHMlo memory B",
                 "IGHMhi_memory_B"="IGHMhi memory B",
                 "pDC"="pDC",
                 "atypical_B"="Atypical B",
                 "CD4\\+_T" = "CD4+ T") # 22 cell types

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC",
                           "#053061"))(300) %>% rev()
#leafcutter
data<-read.delim("/Users/tongyihan/project/AIDA/data/figure/figure1/data/leafcutter_corr_spearson_22.txt",sep = "\t",check.names = F)#spearman
data <- data[which(rownames(data)!="CD4+_T"),which(colnames(data)!="CD4+_T")] # remove CD4+ T 2023.10.27 
corr_leaf <- cor(data[,1:21], use = 'pairwise', method = 'spearman')

#' replace old names from the plot file

for (i in seq_along(label_dict)) {
  colnames(corr_leaf) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],colnames(corr_leaf))
  rownames(corr_leaf) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],rownames(corr_leaf))
}
name_order<-c("pDC","cDC2","CD14+ Monocyte","CD16+ Monocyte","Atypical B","IGHMlo memory B","IGHMhi memory B","Naive B","em CD4+ T",
              "Naive CD4+ T","Naive CD8+ T","cm CD4+ T","Treg","GZMKhi CD8+ T","MAIT","GZMKhi gdT","cyt CD4+ T","GZMBhi CD8+ T","CD56+ NK",
              "CD16+ NK","GZMBhi gdT")
corr_leaf<-corr_leaf[name_order,name_order]

pdf("/Users/tongyihan/project/AIDA/data/figure/figure1/v2/Fig1.leafcutter_spearman_heatmap_without_CD4T_20231027_order_manually.pdf",width = 8,height = 8) # plot figure 2A
#corrplot(corr_leaf,method="color", type='upper', outline = FALSE,tl.srt = 45,tl.col = "black", tl.offset = 0.9, tl.cex = 0.9, cl.pos = 'r',number.font = 6, hclust.method = "complete",order="hclust",
#        col = col2,rect.col="white",rect.lwd=2,addrect=T,addgrid.col="white") #by cluster
corrplot(corr_leaf,method="color", type='upper', outline = FALSE,tl.srt = 45,tl.col = "black", tl.offset = 0.9, tl.cex = 0.9, cl.pos = 'r',number.font = 6, hclust.method = "complete",order="original",
         col = col2,rect.col="white",rect.lwd=2,addrect=T,addgrid.col="white") # order manually
dev.off()


#spliz
mat_spliz<-read.delim("/Users/tongyihan/project/AIDA/data/figure/figure2/new version/summary.celltype_mean.top_variable.1000.figure2A.tsv",check.names = F)
mat_spliz<-column_to_rownames(mat_spliz,var = "gene")
mat_spliz<-mat_spliz[,which(colnames(mat_spliz) %in% cell_type)]
col3 <- colorRampPalette(c("#67001F", "#B2182B",  "#F4A582",
                           "#FFFFFF",  "#92C5DE", "#2166AC",
                           "#053061"))(100) %>% rev()

for (i in seq_along(label_dict)) {
  colnames(mat_spliz) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],colnames(mat_spliz))
}


#name_order<-c("pDC","cDC2","CD14+ Monocyte","CD16+ Monocyte","Atypical B","Naive B","IGHMhi memory B","IGHMlo memory B",
#              "Naive CD4+ T","Naive CD8+ T","em CD4+ T","Treg","cm CD4+ T","GZMKhi gdT","GZMKhi CD8+ T","MAIT","CD56+ NK",
#             "cyt CD4+ T","GZMBhi CD8+ T","CD16+ NK","GZMBhi gdT") # remove CD4+ T,previous one

#name_order<-c("pDC","cDC2","CD14+ Monocyte","CD16+ Monocyte","Atypical B","Naive B","IGHMhi memory B","IGHMlo memory B","em CD4+ T","Treg","cm CD4+ T","Naive CD4+ T","Naive CD8+ T","GZMKhi CD8+ T","GZMKhi gdT",
#"MAIT","CD56+ NK","CD16+ NK","GZMBhi gdT","cyt CD4+ T","GZMBhi CD8+ T") # remove CD4+ T, order is consist with leafcutter
name_order<-c("pDC","cDC2","CD14+ Monocyte","CD16+ Monocyte","Atypical B","IGHMlo memory B","IGHMhi memory B","Naive B","em CD4+ T",
              "Naive CD4+ T","Naive CD8+ T","cm CD4+ T","Treg","GZMKhi CD8+ T","MAIT","GZMKhi gdT","cyt CD4+ T","GZMBhi CD8+ T","CD56+ NK",
              "CD16+ NK","GZMBhi gdT") # consistent with leafcutter
mat_spliz<-mat_spliz[,name_order]
mat_spliz.cor.spearman <- cor(mat_spliz[,c(1:21)], use = 'pairwise', method = 'spearman') %>% as.data.frame()
mat_spliz.cor.spearman<-as.matrix(mat_spliz.cor.spearman)

pdf("/Users/tongyihan/project/AIDA/data/figure/figure2/new version/Fig2A.spliz_spearman_heatmap_without_CD4T20231027_manual_order.pdf",width = 8,height = 8)
corrplot(mat_spliz.cor.spearman,method="color", type='lower', outline = FALSE,tl.srt = 45,tl.col = "black", tl.offset = 0.9, tl.cex = 0.9, cl.pos = 'r',number.font = 6, hclust.method = "complete",order="original",
         col = col3,rect.col="white",rect.lwd=2,addrect=F,addgrid.col="white",col.lim = c(-0.3,1),diag = F,is.corr = TRUE) #by cluster
dev.off()

