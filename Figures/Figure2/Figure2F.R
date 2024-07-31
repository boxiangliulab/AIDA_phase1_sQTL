library(dplyr)
library(ggplot2)
library(ggpubr)

clu2gene <- read.delim('data/aida_clu_gene.txt')
celltypes <- read.table('data/cell_type_use.txt')$V1

flist <- list.files('~/00_tmp/DSG_leafcutter_country_final/diff_leafcutter/')
comps <- unique(gsub('.*country\\.(.*)\\..*\\.txt', '\\1', flist))
DSGs <- data.frame()
tb <- data.frame()
tb_gsea <- list()
res_gsea <- list()
res_df.topIntron.merge <- data.frame()
for (comp in comps) {
  for (celltype in celltypes) {
    res_df1 <- read.delim(paste0('~/00_tmp/DSG_leafcutter_country_final/diff_leafcutter/leafcutter_ds.country.', comp, '.', celltype, '_cluster_significance.txt'))
    res_df2 <- read.delim(paste0('~/00_tmp/DSG_leafcutter_country_final/diff_leafcutter/leafcutter_ds.country.', comp, '.', celltype, '_effect_sizes.txt'))
    
    res_df2$cluster <- gsub(':[0-9]*:[0-9]*:', ':', res_df2$intron)
    
    nrow(res_df2)
    #nrow(res_df2.sum)
    #/ebs1/users/tianchi/project/proj_aida_sqtl/scripts/batch_v1_0/30_DS_pseudobulk_leafcutter/get_filtered_introns.R
    introns_filter <- read.table(paste0('../DSG_leafcutter_filtered_introns/filtered_introns.', celltype, '.txt'))$V1
    
    res_df <- merge(res_df2, res_df1, by = 'cluster', all.x = T)
    res_df <- merge(res_df, clu2gene, by.x = 'cluster', by.y = 'clu')
    res_df <- subset(res_df, intron %in% introns_filter)
    res_df.sum <- res_df %>%
      group_by(cluster) %>%
      summarise(mean_deltapsi = mean(abs(deltapsi)), nrow = n())
    ggplot(res_df.sum, aes(nrow)) + geom_histogram(binwidth = 1, color = 'black') + theme_bw()
    
    # equal delta psi is common
    res_df.topIntron <- res_df %>%
      group_by(cluster) %>%
      top_n(1, abs(deltapsi))

    # significant 
    res_df.topIntron.sig <- subset(res_df.topIntron, p.adjust < .05 & abs(deltapsi) > .1)
    dim(res_df.topIntron.sig)
    res_df.topIntron.sig$celltype <- celltype
    res_df.topIntron.sig$country <- comp
    DSGs <- rbind(DSGs, res_df.topIntron.sig)
    
    res_df.topIntron$celltype <- celltype
    res_df.topIntron$country <- comp
    res_df.topIntron.merge <- rbind(res_df.topIntron.merge, res_df.topIntron)
    
    
    length(unique(res_df.topIntron$genes))
    length(unique(res_df.topIntron.sig$genes))
  }
}
DSGs <- DSGs[order(abs(DSGs$deltapsi), decreasing = T),]
write.table(DSGs, 'data/table.DSGs.country.tsv', quote = F, sep = '\t', row.names = F)
write.table(res_df.topIntron.merge, 'data/table.topIntron.country.tsv', quote = F, sep = '\t', row.names = F)
# post analysis starts from here  
DSGs <- read.delim('data/table.DSGs.country.tsv')
res_df.topIntron.merge <- read.delim('data/table.topIntron.country.tsv')




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

comps_use <- c('Chinese_vs_Indian', 'Chinese_vs_Malay', 'Indian_vs_Malay')
dat.sig <- subset(DSGs, country %in% comps_use & celltype != 'CD4+_T')
dat.sig$chr <- gsub('(.*):.*', '\\1', dat.sig$cluster)
source('data/color_palette.R')
for (i in seq_along(label_dict)) {
  dat.sig$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat.sig$celltype)
  names(celltype_pal) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(celltype_pal))
}
##### 20230826 #####
dat.sig.genes <- dat.sig %>%
  group_by(celltype, genes, country) %>% summarise(count = n())
dat.sig.genes$celltype <- reorder(dat.sig.genes$celltype, dat.sig.genes$celltype, FUN = length)


dat.sig$celltype <- reorder(dat.sig$celltype, dat.sig$celltype, FUN = length)
f2 <- ggplot(dat.sig, aes(x = celltype)) + geom_bar(stat = 'count', position = 'stack', color = 'black', fill = 'white') +
  #   ylab('Count') + ylim(0, round(1.1*max_ct)) +
  ylab('') + ylim(0, round(1.25*max(table(dat.sig$celltype)))) +
  geom_text(inherit.aes = F, aes(x = celltype, label = after_stat(count)), stat = "count", hjust = -.2) +
  theme_classic() + coord_flip() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title.y = element_blank())

country_pal <- c('#61B1B0', '#C7C781', '#93AA9F')
names(country_pal) <- c('Chinese_vs_Malay', 'Chinese_vs_Indian', 'Indian_vs_Malay')
f3 <- ggplot(dat.sig, aes(x = celltype, fill = country)) + geom_bar(stat = 'count', position = 'fill', color = 'black') +
  xlab('Cell type') + ylab('') + geom_abline(slope = 0, intercept = .5, linetype = 2) +
  theme_classic() + scale_fill_manual(values = country_pal) + coord_flip() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y=element_text(size=11,face="bold"))


ggarrange(f3, f2, nrow = 1, align = 'h', legend = 'bottom', widths = c(2,1), common.legend = T)
ggsave('Fig.2f.pdf', width = 5, height = 5)


