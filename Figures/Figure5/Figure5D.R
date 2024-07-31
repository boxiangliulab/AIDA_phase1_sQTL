library(circlize)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
source("/ebs1/users/jingzhixuan/AIDA_sQTL/cell_type_pretty_names.R")

# The final step to draw circos plot
# There will be some manual adjustment to reach the publication quality after
get_var <- function(row) {
        snp <- as.character(row["variant"])
        loc <- strsplit(snp, "_")[[1]][2]
        return(as.integer(loc))
}

get_chr <- function(row) {
        snp <- as.character(row["variant"])
        loc <- strsplit(snp, "_")[[1]][1]
        return(loc)
}
cast <- function(row) {
        chr <- as.character(row["chr"])
        chr <- as.integer(gsub("chr", "", chr))
        return(chr)
}
rou <- function(row) {
        len <- nchar(row["celltype"])
        rou <- anc + (7 - len) * char_len
        return(rou)
}


cells <- c('atypical_B','CD4+_T_cm','CD4+_T_cyt','CD4+_T_em','CD4+_T_naive','CD4+_T','CD8+_T_GZMB+','CD8+_T_GZMK+','CD8+_T_naive','CD14+_Monocyte','CD16+_Monocyte','CD16+_NK','CD56+_NK','cDC2','gdT_GZMK-','gdT_GZMK+','IGHMhi_memory_B','IGHMlo_memory_B','MAIT','naive_B','pDC','Treg') 

circos_path <- "/ebs1/users/jingzhixuan/AIDA_sQTL/summary_overlap.txt"
outdir <- "/ebs1/users/jingzhixuan/AIDA_sQTL/circos_coloc/circos_result/"
trans_list <- "/ebs1/users/jingzhixuan/AIDA_sQTL/all_cells_trans-sQTL/"
dat <- fread(circos_path)
filter <- fread("/ebs1/users/jingzhixuan/AIDA_sQTL/circos_coloc/h4_filtered.txt")
colnames(filter) <- c("gene.y", "gene.x", "celltype", "h4")
filter <- filter[(filter$h4 >= 0.75), ]
ref <- fread("/ebs1/shared/data/reference/hg38/gencode.v32.primary_assembly.annotation.gtf.genelist.tsv")
colnames(ref) <- c("chr", "start1", "start2", "check", "SN", "type", "gene")
ref$pos <- ifelse(ref$check == "+", ref$start1, ref$start2)
color <- read.table('/ebs1/shared/data/aida/ref_palette.new.txt', comment.char = '')
char_len <- 0.06
anc <- 0.7  # 7 characters as standard
dat <- dat %>% select(celltype, variant, pair, gene.x, gene.y)

cell <- inner_join(dat, filter, by = c("celltype", "gene.x", "gene.y"))
colnames(color) <- c("celltype", "color")
cell <- inner_join(cell, color, by = c("celltype"))
print(nrow(cell))
cisgene <- cell %>% select(celltype, variant, gene.y)
cisgene$chr <- apply(cisgene, 1, get_chr)
cisgene$start <- apply(cisgene, 1, get_var)
cisgene$end <- cisgene$start + 1
colnames(cisgene) <- c("celltype", "variant", "label", "chr", "start", "end")
cisgene <- cisgene %>% select(celltype, chr, start, end, label)

transgene <- cell %>% select(celltype, pair, gene.x)
transgene$chr <- apply(transgene, 1, function(row) {
        snp <- as.character(row["pair"])
        chr <- strsplit(snp, ":")[[1]][1]
        return(chr)
})
transgene$start <- apply(transgene, 1, function(row) {
        snp <- as.character(row["pair"])
        loc <- strsplit(snp, ":")[[1]][2]
        return(as.integer(loc))
})
transgene$end <- apply(transgene, 1, function(row) {
        snp <- as.character(row["pair"])
        loc <- strsplit(snp, ":")[[1]][3]
        return(as.integer(loc))
})
colnames(transgene) <- c("celltype", "pair", "label", "chr", "start", "end")
transgene <- transgene %>% select(celltype, chr, start, end, label)
for (i in seq_along(label_dict)) {
 color$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],color$celltype)
}
cisgene$rou <- apply(cisgene, 1, rou)
transgene$rou <- apply(transgene, 1, rou)
r_cis <- ifelse(cisgene$label == "HNRNPLL", 0.9, 0.9)
r_trans <- ifelse(transgene$label == "PTPRC", 0.9, 0.9)
cisgene$rou <- r_cis
transgene$rou <- r_trans
png(paste(outdir, 'circos_plot_total_outsider_divided.png', sep=""), width = 800, height = 800, res = 300)
cisgene <- cisgene %>% select(chr, start, end, label, rou)
transgene <- transgene %>% select(chr, start, end, label, rou)
print(cisgene)
print(transgene)
circos.initializeWithIdeogram(plotType = "ideogram")
circos.genomicLink(cisgene, transgene, col = cell$color, rou1 = cisgene$rou, rou2 = transgene$rou)
circos.clear()
dev.off()
