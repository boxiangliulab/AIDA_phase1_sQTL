library(dplyr)
library(ggplot2)
library(data.table)
library(readr)
library(locuscomparer)
library(cowplot)


# Draw locuscompare plot of all isoforms

# prequists:
# merge_path -- data frame containing SNP with corresponding p-value of both cis-eQTl and trans-sQTL
# ld_path -- ld result from locuscompare website
# total -- rsid binding result from SNPNexus

# Output: scatter plot of locuscompare

sub <- function(row) {
        chr <- as.character(row["VID"])
        chr <- substr(chr, 1, nchar(chr) - 2)
        return(chr)
}

add_label <- function(merged, snp){
    merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, '')
    return(merged)
}

merge_path <- "/ebs1/users/jingzhixuan/AIDA_sQTL/HNRNPLL-PTPRC-intronlized-coloc/checkpoints/"
ld_path <- "/ebs1/users/jingzhixuan/AIDA_sQTL/HNRNPLL-PTPRC-coloc/locuscompare/"
l <- list.files(merge_path)
cells <- c('atypical_B','CD4+_T_cm','CD4+_T_cyt','CD4+_T_em','CD4+_T_naive','CD4+_T','CD8+_T_GZMB+','CD8+_T_GZMK+','CD8+_T_naive','CD14+_Monocyte','CD16+_Monocyte','CD16+_NK','CD56+_NK','gdT_GZMK-','gdT_GZMK+','IGHMhi_memory_B','IGHMlo_memory_B','MAIT','naive_B','pDC','Treg') 

total <- fread("/ebs1/users/jingzhixuan/AIDA_sQTL/HNRNPLL-PTPRC-intronlized-coloc/nexus/rsid_ref.txt")
colnames(total) <- c("VID", "rsid", "chr", "pos", "ref", "alt", "minor", "maf", "contig", "contig_pos", "band")
total$VID <- apply(total, 1, sub)
total <- total %>% select(VID, rsid)
 
face <- data.frame(rsid = character(), pos = integer(), pval1 = double(), pval2 = double(), logp1 = double(), logp2 = double())

for (cell in cells){
  if (paste("chr1:198692373:198703298_", cell,"_merged.txt",sep="") %in% l) {
    # Specific isoform used in paper
    merged <- fread(paste(merge_path,"chr1:198692373:198703298_", cell,"_merged.txt",sep=""))
    merged <- inner_join(merged, total, by = c("VID"))
    ld <- fread(paste(ld_path,cell,"/ld.tsv",sep=""))
    merged <- merged %>% select(rsid, pos, epval, spval, logpe, logps)
    colnames(merged) <- c("rsid", "pos", "pval1", "pval2", "logp1", "logp2")
    merged$chr <- "2"
    merged <- merged[(!duplicated(merged[, c("pos")])), ]
    legend <- TRUE
    legend_position <- c('bottomright','topright','topleft')
    snp <- merged[which.min(merged$pval1*merged$pval2), 'rsid']
    snp = get_lead_snp(merged, snp)
    color <- assign_color(merged$rsid, snp, ld)

    shape <- ifelse(merged$rsid == snp, 23, 21)
    names(shape) <- merged$rsid

    size <- ifelse(merged$rsid == snp, 4, 3)
    names(size) <- merged$rsid

    merged <- add_label(merged, snp)

    p <- make_scatterplot(merged, 'HNRNPLL cis-eQTL', 'PTPRC trans-sQTL', color,
                          shape, size, legend, legend_position)

    ggsave(paste("/ebs1/users/jingzhixuan/AIDA_sQTL/HNRNPLL-PTPRC-intronlized-coloc/all_locus/",cell,"_scatter_raster.pdf",sep=""),p, device = "pdf", width = 4, height = 4)
  }
}
