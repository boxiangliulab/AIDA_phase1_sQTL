#Load qvalue package
suppressMessages(library(qvalue))
library(dplyr)
library(stringr)
library(data.table)

annotation_filter <- fread('~/project/proj_aida_sqtl/data/processed/batch_v1_0/trans_sQTL/filter/gencode.v32.primary_assembly.annotation.pc.lncrna.gtf.genelist.map.0.8.tsv')
annotation <- fread('/ebs1/shared/data/reference/hg38/gencode.v32.primary_assembly.annotation.gtf.genelist.tsv')
crossmap <- fread('/ebs1/shared/data/reference/hg38/mappability/hg38_cross_mappability_strength.txt.gz')

annotation$geneid <- gsub('\\..*', '', annotation$V5)
crossmap$V1 <- gsub('\\..*', '', crossmap$V1)
crossmap$V2 <- gsub('\\..*', '', crossmap$V2)

filter_crossmap <- function(x){
  # chr pos ref alt
  variant <- str_split_fixed(x[1], ':', n = 4)
  target <- x[2]
  target_id <- annotation[annotation$V7 == target, ]$geneid
  nearby_genes <- c()
  for (i in 1:nrow(annotation)) {
    if (annotation[i, 1] == variant[1]){
      if (abs(annotation[i, 2] - as.numeric(variant[2])) < 1000000 | abs(annotation[i, 3] - as.numeric(variant[2])) < 1000000){
        nearby_genes <- c(nearby_genes, annotation[i, ]$geneid)
      }
    }
  }
  crossmap_geneids <- c()
  crossmap_geneids_target <- c()
  for (n in nearby_genes) {
    # only the nearby-target pair
    tmp <- subset(crossmap, V1 == n)
    tmp <- subset(tmp, V2 %in% target_id)
    crossmap_geneids <- unique(c(crossmap_geneids, tmp$V1))
    crossmap_geneids_target <- unique(c(crossmap_geneids_target, tmp$V2))
  }
  crossmap_genelist <- subset(annotation, geneid %in% crossmap_geneids)$V7

  # filter
  if (! target %in% annotation_filter$V7){return('notPCorLNC')}
  else if (length(crossmap_genelist)>0){return(paste0('crossmap','(', paste(crossmap_genelist, collapse = ','), ')'))}
  else{return('PASS')}
}



#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript runFDR_trans_gene.R adjusted.best.txt FDR output.gene.txt"))
cat("\nProcessing QTLtools approximate trans output\n");
cat("  * File best  = [", args[1], "]\n");
cat("  * FDR        = [", args[2], "]\n");
cat("  * Output     = [", args[3], "]\n");


dat <- read.table(args[1], head=FALSE, stringsAsFactors=FALSE)
FDR = as.numeric(args[2])

intron_to_gene <- read.delim('/ebs1/zhangyuntian/aida_clu_gene.txt')
dat$clu <- paste0(lapply(strsplit(dat$V1, ':'), `[[`, 1), ':', lapply(strsplit(dat$V1, ':'), `[[`, 4))

dat.merge <- merge(dat, intron_to_gene)
dat.summary <- dat.merge %>%
 group_by(genes) %>%
 mutate(n = n(), smallest_p = min(V2)) %>%
 slice_min(order_by = V2)
 
# adjust per gene across k phenotypes 
dat.summary$p_gene <- 1-(1-dat.summary$smallest_p)^dat.summary$n

# adjust for multiple phenotypes
dat.summary$qval <- qvalue(dat.summary$p_gene)$qval


##### the last columns are gene, k (# of phenotypes), minimum p-value, p-value per gene, and qvalue #####
write.table(dat.summary, args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
dat.summary.filter <- subset(dat.summary, qval <= FDR)
# filter based on mappability, cross-mappability
dat.summary.filter$mappability <- apply(dat.summary.filter[, c('V4', 'genes')], 1, filter_crossmap)
cat("  * " , nrow(dat.summary.filter) , " are significante out of ", nrow(dat.summary), "\n")
write.table(dat.summary.filter, paste0(args[3], '.filtered.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
