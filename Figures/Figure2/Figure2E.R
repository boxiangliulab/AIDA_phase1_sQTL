library(reshape2)
library(ggplot2)
library(dplyr)

stat <- read.delim('data/table.topIntron.sex.tsv')
rank <- stat %>% subset(genes == 'FLNA') %>% group_by(celltype) %>% slice_min(p.adjust) %>% arrange(-p.adjust) 
ct_list.FLNA <- read.table('data/cell_type_use.QTL.txt')$V1

start.plot <- 154348000
end.plot <- 154356000

df <- data.frame()
celltypes.rank <- c(ct_list.FLNA[!ct_list.FLNA %in% rank$celltype], unique(rank$celltype))
celltypes.use <- c('IGHMlo_memory_B', 'CD8+_T_GZMB+', 'CD4+_T_cyt', 'gdT_GZMK-')
for (ct in celltypes.use) {
  dat.f <- read.table(paste0('data/DS/mergeAll.', ct, '.Female.FLNA_zoom.sampled5000.bam.bdg'))
  dat.f$V2 <- dat.f$V2+1
  dat.f <- subset(dat.f, V2 > start.plot & V3 < end.plot)
  dat.f$celltype <- ct
  
  dat.m <- read.table(paste0('data/DS/mergeAll.', ct, '.Male.FLNA_zoom.sampled5000.bam.bdg'))
  dat.m$V2 <- dat.m$V2+1
  dat.m <- subset(dat.m, V2 > start.plot & V3 < end.plot)
  dat.m$celltype <- ct
  
  dat.f.melt <- melt(dat.f, measure.vars = c(2,3), id.vars = c(1,4,5), value.name = 'checkpoint')
  dat.m.melt <- melt(dat.m, measure.vars = c(2,3), id.vars = c(1,4,5), value.name = 'checkpoint')
  dat.f.melt$sex <- 'Female'
  dat.m.melt$sex <- 'Male'
  
  dat.melt <- rbind(dat.f.melt, dat.m.melt)
  df <- rbind(df, dat.melt)
}

source('src/cell_type_pretty_names.R')
names(label_dict)[which(names(label_dict) == 'naive_B')] <- 'Naive_B'
rank.subset <- subset(rank, celltype %in% celltypes.use)
for (i in seq_along(label_dict)) {
  df$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],df$celltype)
  celltypes.use <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],celltypes.use)
  rank.subset$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],rank.subset$celltype )
}
df$celltype <- factor(df$celltype, levels = celltypes.use)
rank.subset$celltype <- factor(rank.subset$celltype, levels = celltypes.use)
f <- ggplot(df) + geom_line(aes(checkpoint, V4, color = sex)) + theme_classic() +
  scale_color_manual(values = c('#b7352d', '#3c70a4')) + facet_wrap(~celltype, ncol = 1, scales = 'free_x', strip.position = "left") +
  xlim(start.plot, end.plot) +  theme(axis.title = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    strip.text.y.left = element_text(angle=0),strip.background = element_blank(),
                                    strip.placement = "outside") +
  expand_limits(y = 0) + scale_y_continuous(expand = c(0, 0)) +
  geom_text(aes(x = 154354000, y = 800, label = paste0('p.adjust = ', signif(p.adjust, 2))), data = rank.subset)







library(gggenes)
library(data.table)
library(dplyr)
library(ggsci)

gene <- fread("data/DS/FLNA_model.tsv")
features <- data.frame(transcript = 'ENST00000498491',
                       name = 'tss',
                       type = 'tss', 
                       position = 154351203,
                       forward = FALSE)

gene <- subset(gene, start > start.plot & end < end.plot)
gene <- subset(gene, transcript %in% c('ENST00000610817', 
                                       'ENST00000498491',
                                       'ENST00000415241',
                                       'ENST00000422373'))
gene <- gene %>% mutate(exonn = end == 154351203)

p <- ggplot(gene, aes(xmin = start, xmax = end, y = transcript, fill = exonn, color = exonn)) +
  geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"), arrowhead_height = grid::unit(3, "mm")) +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  theme_genes() +
  geom_feature(
    data = features,
    aes(x = position, y = transcript, forward = forward)
  ) +
  xlab("chrX") +
  ylab(NULL) +
  theme(legend.position = "none") +
  xlim(start.plot, end.plot)





library(patchwork)
f + p + plot_layout(ncol = 1, heights = c(5,2))
ggsave('Fig.2e.pdf', width = 7, height = 4)

