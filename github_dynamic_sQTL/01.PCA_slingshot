#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Run PCA and Slingshot on B cells 
# author: Yihan Tong
# date: 2023-01-09
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# conda activate r-4.2.3

#   ____________________________________________________________________________
#   Get option                             ####

library(getopt)
command=matrix(c("file","i",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if( !is.null(args$help) || is.null(args$file))
{ cat(paste(getopt(command,usage=T),"\n"))
  q() }
subdata=args$file
subdata <- readRDS(subdata)
print("check 1 finish reading files")

# .    Import libraries    # 
#pri
library("tidyverse")
#library("dsLib")

#sec
library("Seurat")
library("future")
library("dittoSeq")
library(dplyr)
library(slingshot)
#library(randomcoloR)

# Normalization
DefaultAssay(subdata) <- "integrated"
all.genes <- rownames(subdata)
subdata <- ScaleData(subdata, features = all.genes)
subdata <- RunPCA(subdata, features = VariableFeatures(object = subdata))
### plot seurat clusters
var="seurat_clusters"
pdf("01.pca_reduction_seurat_cluster.pdf",width = 5,height =4)
dittoDimPlot(subdata,var = var,reduction.use = "pca",size = 0.7, main = "",xlab = "PCA_1",
             ylab ="PCA_2", 
             dim.1 = 1,
             dim.2 = 2)
dev.off()

### plot annotation
var="cell_annotation_level3"
pdf("02.pca_reduction_annotation.pdf",width = 5,height = 4)
dittoDimPlot(subdata,var = var,reduction.use = "pca",size = 0.7, main = "",xlab = "PCA_1",
             ylab = "PCA_2", 
             dim.1 = 1,
             dim.2 = 2)
dev.off()

print("check 2 finish plot reduction")
### start pseudotime
em <- Embeddings(subdata, reduction = "pca")
sds <- slingshot(em, clusterLabels = subdata$cell_annotation_level3,start.clus="naive B")
                 #start.clus = "naive B", end.clus = "memory B")
crv1 <- getCurves(sds)
cl<-subdata@meta.data$cell_annotation_level3
cl<-gsub("memory B","#a6cee3",cl)
cl<-gsub("naive B","#fdbf6f",cl)
pdf("03.slingshot_pseudotime_pca.pdf",width=5,height=4)
plot(sds$reducedDim, col = cl, asp = 1, pch = 16,cex=.3)
lines(SlingshotDataSet(crv1), col="black",lwd = 3)
legend('topleft',c("naive B","memory B"),col=c("#fdbf6f","#a6cee3"),title="B cell subtypes")
dev.off()
print("check 3 finish plot curves")

#### plot only one curve
i <- 1
curve_1 <- slingCurves(sds)[[i]]
pt <- curve_1$lambda %>% as.data.frame() %>% set_names("pt")
subdata <- AddMetaData(subdata, pt)
curve_1 <- curve_1$s[curve_1$ord, 1:2]
colnames(curve_1) <- c("PCA_1", "PCA_2")

dittoDimPlot(subdata,var = "pt",reduction.use = "pca",size=0.7,
             dim.1 = 1,
             dim.2 = 2) +
  geom_path(aes(PCA_1, PCA_2), data = as.data.frame(curve_1),  size = 0.7) +
  xlab("PCA 1") +
  ylab("PCA 2") +
  scale_color_viridis_c(name = "Pseudotime")
ggsave("04.curve1.pdf", height = 4, width = 5)

saveRDS(sds,"sds.rds")


pt <- slingPseudotime(sds, na = FALSE) %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode") %>% 
  select(barcode, curve = Lineage1) %>% 
  as_tibble()

md <- subdata[[]] %>% 
  select(barcode) %>% 
  as_tibble()

pt <- full_join(md, pt, by = "barcode")

pt$Q4 <- Hmisc::cut2(pt$curve, 
                     g = 4) %>% 
  `levels<-`(paste0("Q",seq_len(4)))

pt$Q5 <- Hmisc::cut2(pt$curve, 
                     g = 5) %>% 
  `levels<-`(paste0("Q",seq_len(5)))

pt$Q6 <- Hmisc::cut2(pt$curve, 
                     g = 6) %>% 
  `levels<-`(paste0("Q",seq_len(6)))

pt$Q7 <- Hmisc::cut2(pt$curve, 
                     g = 7) %>% 
  `levels<-`(paste0("Q",seq_len(7)))

pt$Q8 <- Hmisc::cut2(pt$curve,
                     g = 8) %>%
  `levels<-`(paste0("Q",seq_len(8)))

pt$Q9 <- Hmisc::cut2(pt$curve,
                     g = 9) %>%
  `levels<-`(paste0("Q",seq_len(9)))

pt$Q10 <- Hmisc::cut2(pt$curve,
                     g = 10) %>%
  `levels<-`(paste0("Q",seq_len(10)))
bins <- pt %>% 
  as.data.frame() %>% 
  column_to_rownames("barcode") %>% 
  select(Q4, Q5, Q6,Q7,Q8,Q9,Q10)
print("check 5 finish quantile")
subdata <- AddMetaData(subdata, bins)
print("check 6 finish add metadata")
quantile<-data.frame(barcode=subdata[[]]$barcode,country=subdata[[]]$country,Q5=subdata[[]]$Q5,Q6=subdata[[]]$Q6,Q7=subdata[[]]$Q7,Q8=subdata[[]]$Q8,Q9=subdata[[]]$Q9,Q10=subdata[[]]$Q10)
write.csv(quantile,"barcode_quantile.csv")
#metadata <- subdata@meta.data
#metadata<-dplyr::filter(metadata,ethnicity=="Japanese")
#write.table(metadata,"metadata_Japanese_only.txt")
saveRDS(subdata,"add_curves_annotation_B_cell_umap.rds")
pdf("05.pca_Q5_quantile.pdf",width = 5,height = 4)
dittoDimPlot(subdata, var = "Q5", reduction.use = "pca",size=0.7,dim.1 = 1,dim.2 = 2)
dev.off()
pdf("05.pca_Q6_quantile.pdf",width = 5,height = 4)
dittoDimPlot(subdata, var = "Q6", reduction.use = "pca",size=0.7,dim.1 = 1,dim.2 = 2)
dev.off()
pdf("05.pca_Q7_quantile.pdf",width = 5,height = 4)
dittoDimPlot(subdata, var = "Q7", reduction.use = "pca",size=0.7,dim.1 = 1,dim.2 = 2)
dev.off()
pdf("05.pca_Q8_quantile.pdf",width = 5,height = 4)
dittoDimPlot(subdata, var = "Q8", reduction.use = "pca",size=0.7,dim.1 = 1,dim.2 = 2)
dev.off()
pdf("05.pca_Q9_quantile.pdf",width = 5,height = 4)
dittoDimPlot(subdata, var = "Q9", reduction.use = "pca",size=0.7,dim.1 = 1,dim.2 = 2)
dev.off()
pdf("05.pca_Q10_quantile.pdf",width = 5,height = 4)
dittoDimPlot(subdata, var = "Q10", reduction.use = "pca",size=0.7,dim.1 = 1,dim.2 = 2)
dev.off()