source('data/cell_type_pretty_names.R')
library(ComplexHeatmap)
celltypes <- read.table('data/cell_type_use.QTL.txt')$V1
color <- read.table('data/ref_palette.new.txt', comment.char = '')
palatte <- color$V2
names(palatte) <- color$V1

dat <- read.delim('data/stats.trans_sQTL.tsv')
for (i in seq_along(label_dict)) {
 dat$cell_type <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat$cell_type)
 celltypes <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],celltypes)
 names(palatte) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(palatte))
}




# using sGenes 
# inferred from later 
setOrder.manual <- c('CD16+ Monocyte', 'CD14+ Monocyte', 'cDC2', 
                     'GZMKhi CD8+ T', 'GZMBhi CD8+ T', 'Naive CD8+ T', 'MAIT',
                     'em CD4+ T', 'cm CD4+ T', 'Naive CD4+ T', 'cyt CD4+ T', 'Treg', 
                     'Naive B', 'IGHMlo memory B', 'IGHMhi memory B', 
                     'CD16+ NK', 'CD56+ NK',
                     'GZMKhi gdT', 'GZMBhi gdT')


lt = list()
#for (ct in celltypes) {
for (ct in setOrder.manual) {
  tmp <- subset(dat, cell_type %in% ct)$gene_id
  lt <- append(lt, list(tmp))
}

names(lt) <- setOrder.manual

m1 = make_comb_mat(lt)
m <- m1[comb_size(m1) >= 3]


# final plot 
library(stringr)
ss = set_size(m)
cs = comb_size(m)


pdf('Fig.5a.pdf', width = 6, height = 5)
UpSet(m, 
      set_order = match(setOrder.manual, set_name(m)), 
      right_annotation = upset_right_annotation(m, gp = gpar(col = NA, fill = palatte[names(ss)])),
      comb_order = order(comb_degree(m), -cs), 
      top_annotation = upset_top_annotation(m, gp = gpar(col = NA, fill = c(rep('black', length(grep(1, comb_degree(m), invert = T))), palatte[set_name(m)]))),
      # only define dot color now (since up and right panel has been changed)
      comb_col = '#00004d'
)
dev.off()


