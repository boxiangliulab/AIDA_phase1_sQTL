library(ggplot2)
library(data.table)
library(dplyr)
source('src/cell_type_pretty_names.R')
source('src/color_palette.R')
celltypes <- read.table('data/cell_type_use.txt')$V1
dat <- read.delim('data/aida_perind.counts.qqnorm.AS_marker.Lynch2004.tsv')

dat.melt <- melt(dat, id.vars = c(1:5, ncol(dat)), value.name = 'normalized_splicing')
# exclude irrelavant intron clusters
unique(dat.melt$cluster)
unique(dat.melt$genes)
dat.melt <- subset(dat.melt, !cluster %in% c('chr1:clu_3479_+', # PTPRC
                                             'chr11:clu_6559_+', # CD44
                                             'chr6:clu_30122_-','chr6:clu_30123_-', # FYN
                                             'chr5:clu_28645_+','chr5:clu_28646_+', # IL7R
                                             'chr16:clu_13772_+','chr16:clu_13773_+', # IL4R
                                             'chr17:clu_15083_-', # PECAM1
                                             'chr10:clu_4311_+' # FAS
                                             ))
# dat.melt <- subset(dat.melt, genes != 'MALT1')

dat.melt$celltype <- gsub('\\.bam', '', substr(dat.melt$variable, 13,100))
dat.melt.group <- dat.melt %>%
  group_by(ID, celltype, genes) %>%
  summarise(average_normalized_splicing = mean(normalized_splicing))


dat.melt.group$celltype <- gsub('\\.\\.1', '-', dat.melt.group$celltype)
dat.melt.group$celltype <- gsub('\\.', '+', dat.melt.group$celltype)
dat.melt.group$iid <- paste(dat.melt.group$genes, dat.melt.group$ID, sep = '__')

for (i in seq_along(label_dict)) {
  dat.melt.group$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat.melt.group$celltype)
  names(celltype_pal) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(celltype_pal))
}

library(reshape2)
dat.melt.group.mtx <- acast(dat.melt.group, iid~celltype, value.var = 'average_normalized_splicing')


ord_level <- c()
for (i in unique(dat.melt.group$genes)) {
  dat.melt.group.mtx.tmp <- acast(subset(dat.melt.group, genes == i), iid~celltype, value.var = 'average_normalized_splicing')
  hclust.tmp <- hclust(dist(dat.melt.group.mtx.tmp, method = "euclidean"), method = "ward.D")
  ord_level <- c(ord_level, rownames(dat.melt.group.mtx.tmp)[hclust.tmp$order])
}
hclust2 <- hclust(dist(t(dat.melt.group.mtx), method = "euclidean"), method = "ward.D")
ord2 <- hclust2$order




# CD44 
g <- 'CD44'
dat.melt.group.tmp <- subset(dat.melt.group, genes == g & celltype %in% grep(' CD4|CD8|MAIT|Treg|gdT', unique(dat.melt.group$celltype), value = T))

dat.melt.group.tmp.mtx <- acast(dat.melt.group.tmp, iid~celltype, value.var = 'average_normalized_splicing')
hclust.tmp <- hclust(dist(dat.melt.group.tmp.mtx, method = "euclidean"), method = "ward.D")
ord_level <- rownames(dat.melt.group.tmp.mtx)[hclust.tmp$order]

ord_level2.df <- subset(dat.melt.group.tmp, iid == 'CD44__chr11:35190065:35214852:clu_6558_+')
ord_level2 <- ord_level2.df$celltype[order(ord_level2.df$average_normalized_splicing, decreasing = T)]


dat.melt.group.tmp$iid <- factor(dat.melt.group.tmp$iid, levels = ord_level)
dat.melt.group.tmp$celltype <- factor(dat.melt.group.tmp$celltype, levels = ord_level2)

iid_label <- c('CD44s',
               'CD44 exon 5-exon v10',
               'CD44 exon v9-exon v10',
               'CD44 exon v10-exon 6',
               'CD44 exon 5-exon v6')
names(iid_label) <- c('CD44__chr11:35190065:35214852:clu_6558_+',
                      'CD44__chr11:35190065:35211246:clu_6558_+',
                      'CD44__chr11:35210054:35211246:clu_6558_+',
                      'CD44__chr11:35211449:35214852:clu_6558_+',
                      'CD44__chr11:35190065:35204512:clu_6558_+')
dat.melt.group.tmp$iid_label <- factor(iid_label[dat.melt.group.tmp$iid], 
                                       levels = iid_label[levels(dat.melt.group.tmp$iid)])
f2 <- ggplot(dat.melt.group.tmp, aes(x=iid_label, y=celltype, fill=average_normalized_splicing)) + geom_tile() +
  scale_fill_viridis_c() + 
  theme(legend.position = 'top',
        axis.text.y = element_text(colour = celltype_pal[ord_level2]),
        axis.text.x = element_text(colour = ifelse(ord_level %in% cluster.sig, 'red', 'grey'), angle = 45, hjust = 1))

ggsave('Fig.2c.pdf', width = 5, height = 6)

