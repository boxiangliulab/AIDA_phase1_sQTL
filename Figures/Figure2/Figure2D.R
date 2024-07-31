library(dplyr)
library(ggplot2)
library(ggpubr)

source('data/color_palette.R')
clu2gene <- read.delim('data/aida_clu_gene.txt')
celltypes <- read.table('data/cell_type_use.txt')$V1

DSGs <- data.frame()
tb <- data.frame()
tb_gsea <- list()
res_gsea <- list()
res_df.topIntron.merge <- data.frame()
for (celltype in celltypes) {
  res_df1 <- read.delim(paste0('~/00_tmp/DSG_leafcutter_sex_revision/diff_leafcutter/leafcutter_ds.sex.', celltype, '_cluster_significance.txt'))
  res_df2 <- read.delim(paste0('~/00_tmp/DSG_leafcutter_sex_revision/diff_leafcutter/leafcutter_ds.sex.', celltype, '_effect_sizes.txt'))
  
  res_df2$cluster <- gsub(':[0-9]*:[0-9]*:', ':', res_df2$intron)
  
  nrow(res_df2)
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
  DSGs <- rbind(DSGs, res_df.topIntron.sig)
  
  res_df.topIntron$celltype <- celltype
  res_df.topIntron.merge <- rbind(res_df.topIntron.merge, res_df.topIntron)
  
  
  length(unique(res_df.topIntron$genes))
  length(unique(res_df.topIntron.sig$genes))
}
write.table(DSGs, 'data/table.DSGs.sex.tsv', quote = F, sep = '\t', row.names = F)
write.table(res_df.topIntron.merge, 'data/table.topIntron.sex.tsv', quote = F, sep = '\t', row.names = F)
# post analysis starts from here  
DSGs <- read.delim('data/table.DSGs.sex.tsv')
res_df.topIntron.merge <- read.delim('data/table.topIntron.sex.tsv')






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

dat.sig <- subset(DSGs, celltype != 'CD4+_T')
dat.sig$chr <- gsub('(.*):.*', '\\1', dat.sig$cluster)
dat.sig <- dat.sig %>%
  mutate(Sign = case_when(deltapsi > 0 ~ "F", deltapsi < 0 ~ "M")) %>% 
  mutate(Chromosome = case_when(grepl("^chr[0-9]{1,}$", chr) ~ 'Autosomal', chr == 'chrX' ~ 'X-linked', chr == 'chrY' ~ 'Y'))
for (i in seq_along(label_dict)) {
  dat.sig$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat.sig$celltype)
  names(celltype_pal) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(celltype_pal))
}
dat.sig$celltype <- reorder(dat.sig$celltype, dat.sig$celltype, FUN = length)

##### update 20230826 #####
dat.sig.genes <- dat.sig %>%
  group_by(celltype, Chromosome, genes) %>% summarise()
dat.sig.genes$celltype <- reorder(dat.sig.genes$celltype, dat.sig.genes$celltype, FUN = length)

library(stringr)
library(ComplexHeatmap)
lt = list()
for (ct in unique(dat.sig.genes$celltype)) {
  tmp <- subset(dat.sig.genes, celltype %in% ct)$genes
  lt <- append(lt, list(tmp))
}
names(lt) <- unique(dat.sig.genes$celltype)
m1 = make_comb_mat(lt)
m <- m1[comb_size(m1) >= 1]
ss = set_size(m)
cs = comb_size(m)

pdf('Fig.2d.pdf', width = 5, height = 5)
UpSet(m, 
      # set_order = order(names(ss), -ss),
      right_annotation = upset_right_annotation(m, gp = gpar(col = 'black', fill = 'black')),
      comb_order = order(comb_degree(m), -cs),
      top_annotation = upset_top_annotation(m, gp = gpar(col = 'black', fill = 'black')),
      # only define dot color now (since up and right panel has been changed)
      comb_col = '#00004d'
)
dev.off()

