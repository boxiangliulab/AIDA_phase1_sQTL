library(ggplot2)
library(data.table)
library(dplyr)
library(reshape2)

source('src/cell_type_pretty_names.R')
source('src/color_palette.R')
celltypes <- read.table('data/cell_type_use.txt')$V1
dat <- read.delim('data/aida_perind.counts.qqnorm.AS_marker.Lynch2004.ratio.tsv')
dat.melt <- melt(dat, id.vars = c(1, 2, ncol(dat)), value.name = 'normalized_splicing')
# exclude irrelavant intron clusters
dat.melt <- subset(dat.melt, !cluster %in% c('chr1:clu_3479_+', # PTPRC
                                             'chr11:clu_6559_+' # CD44
))
colnames(dat.melt)[2] <- 'ID'

dat.melt$celltype <- gsub('\\.bam', '', substr(dat.melt$variable, 13,100))
dat.melt$normalized_splicing <- sapply(dat.melt[,5], function(x) eval(parse(text = x)))
dat.melt <- na.omit(dat.melt)
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


cluster.sig <- c('PTPRC__chr1:198692373:198703298:clu_3480_+',
                 'CD44__chr11:35190065:35214852:clu_6558_+')


# PTPRC 
g <- 'PTPRC'
dat.melt.group.tmp <- subset(dat.melt.group, iid != 'PTPRC__chr1:198692373:198694018:clu_3480_+' &
                               genes == g & celltype %in% grep(' CD4|CD8|MAIT|Treg|gdT', unique(dat.melt.group$celltype), value = T))

iid_label <- c('PTPRC: exon6-exon7 (RC+)',
               'PTPRC: exon4-exon5 (RAB+)',
               'PTPRC: exon5-exon6 (RBC+)',
               'PTPRC: exon3-exon4 (RA+)',
               'PTPRC: exon3-exon7 (RO)',
               'PTPRC: exon3-exon5 (RB+)',
               'PTPRC: exon5-exon7 (RB+)')
names(iid_label) <- c('PTPRC__chr1:198702530:198703298:clu_3480_+',
                      'PTPRC__chr1:198696909:198699564:clu_3480_+',
                      'PTPRC__chr1:198699704:198702387:clu_3480_+',
                      'PTPRC__chr1:198692373:198696712:clu_3480_+',
                      'PTPRC__chr1:198692373:198703298:clu_3480_+',
                      'PTPRC__chr1:198692373:198699564:clu_3480_+',
                      'PTPRC__chr1:198699704:198703298:clu_3480_+')

dat.melt.group.tmp$iid_label <- iid_label[dat.melt.group.tmp$iid]

dat.melt.group.tmp.mtx <- acast(dat.melt.group.tmp, iid~celltype, value.var = 'average_normalized_splicing')
hclust.tmp <- hclust(dist(dat.melt.group.tmp.mtx, method = "euclidean"), method = "ward.D")
ord_level <- rownames(dat.melt.group.tmp.mtx)[hclust.tmp$order]

ord_level2 <- c('Naive CD8+ T', 
                'Naive CD4+ T',
                'GZMBhi CD8+ T',
                'GZMBhi gdT',
                'GZMKhi CD8+ T',
                'GZMKhi gdT', 
                'MAIT',
                'Treg',
                'cyt CD4+ T',
                'em CD4+ T',
                'cm CD4+ T')

dat.melt.group.tmp$iid <- factor(dat.melt.group.tmp$iid, levels = ord_level)
dat.melt.group.tmp$iid_label <- factor(dat.melt.group.tmp$iid_label, 
                                       levels = iid_label[levels(dat.melt.group.tmp$iid)])
dat.melt.group.tmp$celltype <- factor(dat.melt.group.tmp$celltype, levels = ord_level2)

## 
dat.melt.group.tmp.dcast <- dcast(dat.melt.group.tmp, celltype~iid_label, value.var = 'average_normalized_splicing')
dat.melt.group.tmp.dcast.ratio <- dat.melt.group.tmp.dcast
for (i in 1:nrow(dat.melt.group.tmp.dcast.ratio)) {
  dat.melt.group.tmp.dcast.ratio[i,3] <- log2(dat.melt.group.tmp.dcast.ratio[i,3]/sum(dat.melt.group.tmp.dcast[i, c(2, 4:8)]))
  for (j in c(2, 4:8)) {
    dat.melt.group.tmp.dcast.ratio[i,j] <- log2(dat.melt.group.tmp.dcast.ratio[i,j]/dat.melt.group.tmp.dcast[i,3])
  }
}
dat.melt.group.tmp.dcast.tmp <- melt(dat.melt.group.tmp.dcast.ratio, id.vars = 1, variable.name = 'iid_label', value.name = 'average_normalized_splicing')

# update x cluster
hclust.tmp <- hclust(dist(t(as.matrix(dat.melt.group.tmp.dcast.ratio[,-1])), method = "euclidean"), method = "ward.D")
ord_level <- colnames(dat.melt.group.tmp.dcast)[2:8][hclust.tmp$order]
# manually 
ord_level <- c('PTPRC: exon3-exon7 (RO)',
               'PTPRC: exon6-exon7 (RC+)',
               'PTPRC: exon4-exon5 (RAB+)',
               'PTPRC: exon5-exon6 (RBC+)',
               'PTPRC: exon3-exon4 (RA+)',
               'PTPRC: exon5-exon7 (RB+)',
               'PTPRC: exon3-exon5 (RB+)')
dat.melt.group.tmp.dcast.tmp$iid_label <- factor(dat.melt.group.tmp.dcast.tmp$iid_label, 
                                                 levels = ord_level)

f1 <- ggplot(dat.melt.group.tmp.dcast.tmp, aes(x=iid_label, y=celltype, fill=average_normalized_splicing)) + geom_tile() +
  scale_fill_viridis_c() + 
  theme(legend.position = 'top',
        # axis.text.y = element_text(colour = celltype_pal[colnames(dat.melt.group.tmp.mtx)[ord2]]),
        axis.text.y = element_text(colour = celltype_pal[as.vector(ord_level2)]),
        axis.text.x = element_text(colour = ifelse(ord_level %in% iid_label[cluster.sig], 'red', 'grey'), angle = 45, hjust = 1))
ggsave('Fig.2b.pdf', f1, width = 6, height = 6)



