### This file was used to generate Figure3A of this paper and make a summary to the independent cis-sQTLs we identified previously.

library(tidyverse)
library(data.table)
library(fs)
library(stringr)
library(dplyr)
library(ggpmisc)
source("utils.R")

color_cell<-data.frame(cell_type=c("CD14+_Monocyte","CD16+_Monocyte","CD16+_NK","CD4+_T_cm","CD4+_T_cyt","CD4+_T_em","CD4+_T_naive","CD56+_NK",
                              "CD8+_T_GZMB+","CD8+_T_GZMK+","CD8+_T_naive","IGHMhi_memory_B","IGHMlo_memory_B","MAIT","Treg","cDC2",
                              "gdT_GZMK+","gdT_GZMK-","naive_B"),
                       color=c("#fee090","#fdae61", "#bf812d","#9e9ac8","#6a51a3","#807dba",
                               "#bcbddc","#8c510a","#3690c0","#74a9cf","#a6bddb", "#a1d99b","#74c476",
                               "#0570b0","#CCCC4D","#f46d43","#f768a1", "#dd3497","#c7e9c0"))

donor_list_dir <- "/ebs1/users/tongyihan/AIDA/data/14.individual_cell_num/cell_type_donor_list/full_file/"
out_dir <- "sqtl_summary/"
dir.create(out_dir, showWarnings = FALSE)

# Read gene annotation
gene_annotation <- read_gene_annotation(gene_annotation_path)

# Donor list
donor_list <- read_donor_list(donor_list_dir, cell_types = cell_types)
cell_type_summary <- donor_list %>%
    group_by(cell_type) %>%
    summarize(
        n_donors = n(),
        mean_cells = mean(number),
        median_cells = median(number),
        total_cells = sum(number),
    )

# Read permutation pass QTL
perm_qtl <- read_perm_qtl_all_cell_types(cell_types, data_dir = perm_dir)
perm_qtl_2 <- split_multiple_genes(perm_qtl)
perm_qtl_3 <- merge(perm_qtl_2,
    gene_annotation[, c("gene_name", "gene_id", "type")],
    on = "gene_name"
)

# Read conditional pass QTL
indep_qtl <- read_indep_qtl_all_cell_types(cell_types, conditional_dir)
indep_qtl_2 <- split_multiple_genes(indep_qtl)
indep_qtl_3 <- merge(indep_qtl_2,
    gene_annotation[, c("gene_name", "gene_id", "type")],
    on = "gene_name"
)

# Check duplicates
indep_qtl_3 %>%
    select(gene_name, variant_id, cell_type) %>%
    unique() %>%
    nrow() # 15038
indep_qtl_3 %>% nrow() # 15106
indep_qtl_3[indep_qtl_3 %>%
    select(gene_name, variant_id, cell_type) %>%
    duplicated(), ]
indep_qtl_3 %>% filter(gene_name == "AHRR", cell_type == "CD14+_Monocyte", variant_id == "chr5:303223:A:G")
indep_qtl_3 %>%
    select(-gene_id) %>%
    unique() %>%
    nrow() # 15080
indep_qtl_4 <- indep_qtl_3 %>%
    select(-gene_id) %>%
    unique()
indep_qtl_4[indep_qtl_4 %>%
    select(gene_name, variant_id, cell_type) %>%
    duplicated(), ]
indep_qtl_4 %>% filter(
    gene_name == "GGACT",
    variant_id == "chr13:100588735:C:A",
    cell_type == "CD16+_Monocyte"
) # same gene_name, different intron clusters.
# Checked UCSC to confirm both intron clusters below to GGACT.
# Check number of sQTL by gene types:
indep_qtl_3 %>%
    group_by(type) %>%
    summarize(n = n())
#    type                                   n
#    <chr>                              <int>
#  1 TR_C_gene                             41
#  2 TR_J_gene                            437
#  3 TR_V_gene                              6
#  4 lncRNA                              1224
#  5 misc_RNA                              36
#  6 polymorphic_pseudogene                22
#  7 protein_coding                     12569
#  8 transcribed_processed_pseudogene      84
#  9 transcribed_unitary_pseudogene         2
# 10 transcribed_unprocessed_pseudogene   600
# 11 unprocessed_pseudogene                85

# What proportion of tested genes harbor QTLs?
perm_qtl_3 %>%
    select(gene_name, type) %>%
    unique() %>%
    group_by(type) %>%
    summarize(n = n())
#    type                                   n
#    <chr>                              <int>
#  1 IG_C_gene                              1
#  2 IG_J_gene                              5
#  3 IG_V_gene                              6
#  4 TR_C_gene                              3
#  5 TR_J_gene                             49
#  6 TR_V_gene                              2
#  7 lncRNA                               313
#  8 misc_RNA                               7
#  9 polymorphic_pseudogene                 1
# 10 protein_coding                      6552
# 11 transcribed_processed_pseudogene      11
# 12 transcribed_unitary_pseudogene         6
# 13 transcribed_unprocessed_pseudogene    87
# 14 unprocessed_pseudogene                12

(2581 + 174) / (313 + 6552)
# 0.401311

indep_qtl_3 %>%
    select(gene_name, type) %>%
    unique() %>%
    group_by(type) %>%
    summarize(n = n())
#    type                                   n
#    <chr>                              <int>
#  1 TR_C_gene                              3
#  2 TR_J_gene                             49
#  3 TR_V_gene                              1
#  4 lncRNA                               174
#  5 misc_RNA                               4
#  6 polymorphic_pseudogene                 1
#  7 protein_coding                      2581
#  8 transcribed_processed_pseudogene       8
#  9 transcribed_unitary_pseudogene         2
# 10 transcribed_unprocessed_pseudogene    51
# 11 unprocessed_pseudogene                 6

# What proportion of tested genes harbor QTLs?
perm_qtl_3 %>%
    select(gene_name, type) %>%
    unique() %>%
    group_by(type) %>%
    summarize(n = n())
#    type                                   n
#    <chr>                              <int>
#  1 IG_C_gene                              1
#  2 IG_J_gene                              5
#  3 IG_V_gene                              6
#  4 TR_C_gene                              3
#  5 TR_J_gene                             49
#  6 TR_V_gene                              2
#  7 lncRNA                               313
#  8 misc_RNA                               7
#  9 polymorphic_pseudogene                 1
# 10 protein_coding                      6552
# 11 transcribed_processed_pseudogene      11
# 12 transcribed_unitary_pseudogene         6
# 13 transcribed_unprocessed_pseudogene    87
# 14 unprocessed_pseudogene                12

(2581 + 174) / (313 + 6552)
# 0.401311

indep_qtl_3 %>%
    select(gene_name, type) %>%
    unique() %>%
    group_by(type) %>%
    summarize(n = n())
#    type                                   n
#    <chr>                              <int>
#  1 TR_C_gene                              3
#  2 TR_J_gene                             49
#  3 TR_V_gene                              1
#  4 lncRNA                               174
#  5 misc_RNA                               4
#  6 polymorphic_pseudogene                 1
#  7 protein_coding                      2581
#  8 transcribed_processed_pseudogene       8
#  9 transcribed_unitary_pseudogene         2
# 10 transcribed_unprocessed_pseudogene    51
# 11 unprocessed_pseudogene                 6

# Make plot data
plot_data_1 <- indep_qtl_3 %>%
    group_by(cell_type, gene_name) %>%
    summarize(sqtl = n())

plot_data_2 <- plot_data_1 %>%
    mutate(sqtl = ifelse(sqtl >= 5, ">=5", as.character(sqtl))) %>%
    mutate(sqtl = factor(sqtl, levels = c("1", "2", "3", "4", ">=5")))

plot_data_3 <- plot_data_2 %>%
    group_by(cell_type) %>%
    summarize(sgene = n())

cell_type_order <- plot_data_3 %>%
    arrange(desc(sgene)) %>%
    select(cell_type) %>%
    unlist() %>%
    unname()

norm_const <- max(plot_data_3$sgene) * 1.10

custom_sqtl_color <- c(
    "1" = "lightblue", "2" = "lightblue1",
    "3" = "lightblue2", "4" = "lightblue3",
    ">=5" = "lightblue4"
)
newcellname<-c("Naive CD4+ T","cm CD4+ T","CD16+ NK","CD14+ Monocyte","Naive CD8+ T","GZMB+ CD8+ T","CD16+ Monocyte",
           "GZMK+ CD8+ T","em CD4+ T","IGHMlo memory B","Naive B","cDC2","MAIT","Treg",
           "IGHMhi memory B","GZMK+ gdT","cyt CD4+ T","GZMK- gdT","CD56+ NK")
# Make plot: number of independent QTL per cell type?
p1 <- ggplot() +
    geom_bar(
        data = plot_data_2,
        aes(x = cell_type, fill = sqtl),
        position = "fill"
    ) +
    geom_point(
        data = plot_data_3,
        aes(x = cell_type, y = sgene / norm_const),
        color = "red2"
    ) +
    scale_y_continuous(
        name = "Proportion of sGenes",
        expand = c(0.005, 0.005),
        sec.axis = sec_axis(
            trans = ~ . * norm_const,
            name = "Number of sGenes",
            breaks = c(0, 300, 600, 900, 1200, 1500)
        )
    ) +
    scale_x_discrete(limits = cell_type_order) +
    scale_fill_manual(
        values = custom_sqtl_color,
        name = "Independent\nsQTL/sGene",
        labels = c("1", "2", "3", "4", expression("" >= "5"))
    ) +
    xlab(NULL) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(
            fill = alpha("grey90", 0.8),
            color = "black",
            linewidth = 0.5, linetype = "solid"
        ),
        axis.line.y.right = element_line(color = "red2"),
        axis.ticks.y.right = element_line(color = "red2"),
        axis.text.y.right = element_text(color = "red2"),
        axis.title.y.right = element_text(color = "red2")
    )

ggsave(file.path(out_dir, "independent_qtl.pdf"), p1, height = 4, width = 6)
ggsave(file.path(out_dir, "independent_qtl.png"), p1, height = 4, width = 6)

# What proportion of genes have more than 1 QTL?
multi_sqtl_gene <- plot_data_1 %>%
    group_by(cell_type) %>%
    summarize(
        multi_sqtl_sgene = sum(sqtl > 1),
        total_sgene = n(),
        prop_multi_sqtl_gene = multi_sqtl_sgene / total_sgene
    ) %>%
    arrange(desc(prop_multi_sqtl_gene))
multi_sqtl_gene <- merge(multi_sqtl_gene, cell_type_summary, on = "cell_type")
formula <- y ~ x
multi_sqtl_gene<-merge(multi_sqtl_gene,color_cell,by="cell_type")
print(multi_sqtl_gene)

p2 <- ggplot(multi_sqtl_gene, aes(x = total_sgene, y = prop_multi_sqtl_gene)) +
    geom_point(aes(color=cell_type),size=5) +
    scale_color_manual(values=multi_sqtl_gene$color)+
    stat_poly_eq(use_label(c("adj.R2")), formula = formula, label.x = "right", label.y = "bottom") +
    stat_poly_line(formula = formula) +
    theme_classic() +
    xlab("Number of sGenes") +
    ylab("Proportion of sGenes\nwith >1 sSNPs")+guides(color=FALSE)

ggsave(file.path(out_dir, "sgene_vs_allele_heterogeneity.pdf"), height = 3, width = 3)
