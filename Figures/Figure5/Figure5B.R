library(tidyverse)
library(data.table)
library(ggpmisc)
source("src/qtl_summary/utils.R")
source("src/trans_qtl_summary/utils.R")
library(dplyr)

col_names <- c("chr", "start", "end", "strand", "gene_id", "type", "gene_name")
gene_annotation <- fread('data/gencode.v32.primary_assembly.annotation.gtf.genelist.tsv', col.names = col_names)
gene_annotation <- gene_annotation %>% mutate(gene_id = str_split_i(gene_id, fixed("."), 1))

mappability <- fread("data/hg38_gene_mappability.txt.gz")
colnames(mappability) <- c("gene_id", "mappability")
mappability <- mappability %>% mutate(gene_id = str_split_i(gene_id, fixed("."), 1))


trans_dir <- 'data/QTL_trans'
trans <- read_trans(cell_types, trans_dir)
trans <- trans %>% rename(gene_name = gene)
trans_filt_1 <- trans %>% filter(qval < 0.01)
trans_filt_2 <- merge(trans_filt_1, gene_annotation[, c("gene_name", "gene_id", "type")], by = "gene_name")
trans_filt_3 <- merge(trans_filt_2[, -"mappability"], mappability, by = "gene_id")
trans_filt_4 <- trans_filt_3 %>% filter(mappability == 1)

conditional_dir <- 'data/QTL_cis'
indep_qtl <- read_indep_qtl_all_cell_types(cell_types, conditional_dir)
indep_qtl_2 <- split_multiple_genes(indep_qtl)
indep_qtl_3 <- merge(indep_qtl_2,
                     gene_annotation[, c("gene_name", "gene_id", "type")],
                     on = "gene_name"
)

cis_per_cell_type <- indep_qtl_3 %>%
  group_by(cell_type) %>%
  summarize(n_cis = length(unique(gene_name)))

trans_per_cell_type <- trans_filt_3 %>%
  group_by(cell_type) %>%
  summarize(n_trans = length(unique(gene_id)))

trans_sharing <- trans_filt_3 %>%
  group_by(gene_id) %>%
  summarize(n_sGene = n()) %>%
  arrange(desc(n_sGene))


donor_list_dir <- 'data/cell_type_donor_list/full_file'
donor_list <- read_donor_list(donor_list_dir, cell_types = cell_types)
cell_type_summary <- donor_list %>%
  group_by(cell_type) %>%
  summarize(
    n_donors = n()
  )

plot_1 <- merge(cis_per_cell_type, trans_per_cell_type, by = "cell_type")
plot_2 <- merge(plot_1, cell_type_summary)

color <- read.table('data/ref_palette.new.txt', comment.char = '')
colors <- color$V2
names(colors) <- color$V1

f <- ggplot(plot_2, aes(x = n_donors, y = n_trans, group = 1, fill = cell_type, size = n_cis)) +
  geom_point(pch=21) + guides(fill="none") +
  theme_classic() +
  stat_correlation(method = "spearman") +
  xlab("Number of donors") +
  ylab("Number of trans sGenes") +
  scale_fill_manual(name = "Cell types", values = colors)
ggsave('Fig.5b.pdf', f, width=4, height=3.5)
