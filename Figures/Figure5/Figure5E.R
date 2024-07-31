library(ggplot2)
library(ggpubr)
library(patchwork)

source('src/cell_type_pretty_names.R')
source('src/color_palette.R')
celltypes <- read.table('data/cell_type_use.QTL.txt')$V1

expr <- read.delim('data/hnRNPLL.expr.celltype.tsv')
expr2 <- read.delim('data/PTPRC.expr.celltype.tsv')
expr <- subset(expr, celltype %in% celltypes)
expr2 <- subset(expr2, celltype %in% celltypes)

for (i in seq_along(label_dict)) {
  expr$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],expr$celltype)
  expr2$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],expr2$celltype)
  names(celltype_pal) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(celltype_pal))
}

celltype_seq <- unlist(label_dict)
expr$celltype <- factor(expr$celltype, levels = celltype_seq)
expr2$celltype <- factor(expr2$celltype, levels = celltype_seq)


H3 <- read.table('data/h3_all.txt', header = T)
H4 <- read.table('data/h4_all.txt', header = T)
for (i in seq_along(label_dict)) {
  H3$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]], H3$celltype)
  H4$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]], H4$celltype)
}
H4$celltype <- factor(H4$celltype, levels = celltype_seq)


##### v4 #####
cis <- read.delim('data/HNRNPLL.eigenMT.txt')
## trans-sQTL does not necessarily to have the strongest effect on chr1:198692373:198703298:clu_3480_+
trans <- read.delim('data/PTPRC.trans_p.txt')
trans <- subset(trans, celltype != 'IGHMlo_memory_B')

cis$cistrans <- 'cis-eQTL'
cis$logp <- -log10(cis$BF)
trans$cistrans <- 'trans-sQTL'
# trans$logp <- trans$logps
trans$logp <- -log10(trans$p_gene)
for (i in seq_along(label_dict)) {
  cis$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],cis$celltype)
  trans$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],trans$celltype)
  names(celltype_pal) <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],names(celltype_pal))
}
cistrans <- rbind(cis[, c('logp', 'celltype', 'cistrans')], trans[, c('logp', 'celltype', 'cistrans')])

celltype_seq.h4 <- H4$celltype[order(H4$h4)]
celltype_seq.cis <- cis$celltype[order(cis$logp)]
cistrans$celltype <- factor(cistrans$celltype, levels = celltype_seq.cis)
fcistrans_cis <- ggplot(subset(cistrans, cistrans == 'cis-eQTL'), aes(x=celltype, y=cistrans, fill=logp)) + geom_tile() + theme_classic() +
  scale_fill_viridis_c(option = 'plasma') + xlab('') + coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylab(bquote(.(title) ~ -log[10]*'(P)'))
fcistrans_trans <- ggplot(subset(cistrans, cistrans == 'trans-sQTL'), aes(x=celltype, y=cistrans, fill=logp)) + geom_tile() + theme_classic() +
  scale_fill_viridis_c(option = 'plasma') + xlab('') + coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylab(bquote(.(title) ~ -log[10]*'(P)'))


H3$celltype <- factor(H3$celltype, levels = celltype_seq.cis)
H4$celltype <- factor(H4$celltype, levels = celltype_seq.cis)
f3 <- ggplot(H4, aes(x=celltype, y=h4, fill=celltype)) + geom_bar(stat = 'identity') + theme_classic() +
  scale_fill_manual(values = celltype_pal) + xlab('') + geom_hline(yintercept = .75, linetype = 2) +
  theme(axis.ticks.y = element_blank(), legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1),) + coord_flip() +
  scale_y_continuous(breaks = c(0,.75))
f4 <- ggplot(H3, aes(x=celltype, y=h3, fill=celltype)) + geom_bar(stat = 'identity') + theme_classic() +
  scale_fill_manual(values = celltype_pal) + xlab('') + 
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1),) + coord_flip() +
  scale_y_continuous(breaks = c(0,1), limits = c(0,1))

f3 + f4 + fcistrans_cis + fcistrans_trans + plot_layout(nrow = 1)
ggsave('Fig.5e.pdf', width = 5, height = 4)

