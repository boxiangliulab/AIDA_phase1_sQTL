##This script was used to generate Figure3D,E,F and compare the relationship between number of sGenes and number of donors, junction read counts; we also compared  
##the relationship between Proportion of sGenes with >1 sSNPs and number of sGenes.

library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggpmisc)
source("utils.R")

color_cell<-data.frame(cell_type=c("CD14+_Monocyte","CD16+_Monocyte","CD16+_NK","CD4+_T_cm","CD4+_T_cyt","CD4+_T_em","CD4+_T_naive","CD56+_NK",
                              "CD8+_T_GZMB+","CD8+_T_GZMK+","CD8+_T_naive","IGHMhi_memory_B","IGHMlo_memory_B","MAIT","Treg","cDC2",
                              "gdT_GZMK+","gdT_GZMK-","naive_B"),
                       color=c("#fee090","#fdae61", "#bf812d","#9e9ac8","#6a51a3","#807dba",
                               "#bcbddc","#8c510a","#3690c0","#74a9cf","#a6bddb", "#a1d99b","#74c476",
                               "#0570b0","#CCCC4D","#f46d43","#f768a1", "#dd3497","#c7e9c0"))

# Constants
out_dir <- "picture/"
dir.create(out_dir, showWarnings = FALSE)

# Read gene annotation
gene_annotation <- read_gene_annotation(gene_annotation_path)

# Read conditional pass QTL
indep_qtl <- read_indep_qtl_all_cell_types(cell_types, conditional_dir)
indep_qtl_2 <- split_multiple_genes(indep_qtl)
indep_qtl_3 <- merge(indep_qtl_2,
    gene_annotation[, c("gene_name", "gene_id", "type")],
    on = "gene_name"
)
indep_qtl_4 <- indep_qtl_3 %>%
    group_by(cell_type, gene_name) %>%
    summarize(sqtl = n(), fore_slope = max(abs(fore_slope))) # sQTL per sGene
indep_qtl_5 <- indep_qtl_4 %>%
    group_by(cell_type) %>%
    summarize(
        sgene = n(),
        mean_slope = mean(fore_slope),
        median_slope = median(fore_slope)
    ) # sGene per cell type

# Read donor list
donor_list <- read_donor_list(donor_list_dir, cell_types = cell_types)

# Read numerators (intron counts)
numerator <- fread(numerator_file)
colnames(numerator)[1] <- "chrom"
col_names <- colnames(numerator)
col_data <- data.table(col_names = col_names[2:length(col_names)])
col_data <- col_data %>% mutate(
    DCP_ID = str_split_i(col_names, fixed("."), 1),
    cell_type = str_split_i(col_names, fixed("."), 2)
)

# Calculate library size
library_size <- data.frame(colSums(numerator[, -1])) %>%
    as.data.table(keep.rownames = TRUE)
colnames(library_size) <- c("col_names", "library_size")
library_size_2 <- merge(library_size, col_data, on = "col_names")
donor_list_2 <- merge(donor_list, library_size_2, on = c("DCP_ID", "cell_type"))

# Calculate proportion of detected introns (only chr22 for speed)
numerator_chr22 <- numerator %>% filter(startsWith(chrom, "chr22"))
numerator_chr22_long <- numerator_chr22 %>%
    pivot_longer(!chrom, names_to = "col_names", values_to = "counts")
numerator_chr22_long <- merge(numerator_chr22_long, col_data, on = "col_name")
# detected_chr22 <- numerator_chr22_long %>%
#     group_by(cell_type, chrom) %>%
#     summarize(
#         detected = sum(counts > 0),
#         total = n(),
#         prop_detected = detected / total
#     )
# mean_detected_chr22 <- detected_chr22 %>%
#     group_by(cell_type) %>%
#     summarize(
#         mean_detected = mean(detected),
#         mean_prop_detected = mean(prop_detected)
#     )
detected_chr22 <- numerator_chr22_long %>%
    group_by(cell_type, DCP_ID) %>%
    summarize(
        detected = sum(counts > 0),
        total = n(),
        prop_detected = detected / total
    )
mean_detected_chr22 <- detected_chr22 %>%
    group_by(cell_type) %>%
    summarize(
        mean_detected = mean(detected),
        mean_prop_detected = mean(prop_detected)
    )
# Summary per cell type
cell_type_summary <- donor_list_2 %>%
    group_by(cell_type) %>%
    summarize(
        n_donors = n(),
        mean_cells = mean(number),
        median_cells = median(number),
        total_cells = sum(number),
        mean_libsize = mean(library_size),
        median_libsize = median(library_size),
        total_libsize = sum(library_size)
    )

cell_type_summary_1 <- merge(cell_type_summary, indep_qtl_5, on = "cell_type")
cell_type_summary_2 <- merge(cell_type_summary_1, mean_detected_chr22, on = "cell_type")
##add
cell_type_summary_1<-merge(cell_type_summary_1,color_cell,by="cell_type")
# Make plot: number of donors vs number of sGenes
formula_1 <- y ~ x

p1 <- ggplot(cell_type_summary_1, aes(x = n_donors, y = sgene)) +
    geom_point(aes(color=cell_type),size=5) +
    scale_color_manual(values=cell_type_summary_1$color)+
    scale_y_log10() +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_1) +
    stat_poly_line(formula = formula_1) +
    theme_classic() +
    xlab("Number of donors") +
    ylab("Number of sGenes") + guides(color=FALSE)
ggsave(file.path(out_dir, "sgene_vs_donor.pdf"), p1, height = 3, width = 3)


formula_2 <- y ~ x
print(cell_type_summary_1)
p2 <- ggplot(cell_type_summary_1, aes(x = mean_cells, y = sgene)) +
    geom_point() +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_2) +
    stat_poly_line(formula = formula_2) +
    theme_classic() +
    xlab("Cells per donor") +
    ylab("Number of sGenes")
ggsave(file.path(out_dir, "sgene_vs_mean_cells.png"), p2, height = 3, width = 3)


p3 <- ggplot(donor_list_2, aes(x = number, y = library_size)) +
    geom_point() +
    stat_poly_eq(use_label(c("adj.R2")), formula = y ~ x + 0) +
    stat_poly_line(formula = y ~ x + 0) +
    xlab("Number of cells") +
    ylab("Library size") +
    theme_classic()
ggsave(file.path(out_dir, "library_size_vs_number_of_cells.png"),
    p3,
    height = 3, width = 3.5
)
formula_4 <- y ~ x
p4 <- ggplot(cell_type_summary_1, aes(x = mean_libsize, y = sgene)) +
    geom_point(aes(color=cell_type),size=5) +
    scale_color_manual(values=cell_type_summary_1$color)+
    theme_classic() +
    xlab("Junction read counts") +
    ylab("Number of sGenes") +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_4) +
    stat_poly_line(formula = formula_4) + guides(color=FALSE)
ggsave(file.path(out_dir, "sgene_vs_mean_read_count.pdf"),
    p4,
    height = 3, width = 3
)

# Make plot: sample size vs effect size
formula_5 <- y ~ x + poly(x, 2)
p5 <- ggplot(cell_type_summary_1, aes(x = n_donors, y = mean_slope)) +
    geom_point() +
    xlab("Number of donors") +
    ylab("Mean absolute effect size") +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_5, label.x = "right") +
    stat_poly_line(formula = formula_5) +
    theme_classic()
ggsave(file.path(out_dir, "effect_size_vs_sample_size.png"), p5, height = 3, width = 3)

# Make plot: library size vs effect size
formula_6 <- y ~ x
p6 <- ggplot(cell_type_summary_1, aes(x = mean_libsize, y = mean_slope)) +
    geom_point() +
    xlab("Library size") +
    ylab("Mean absolute effect size") +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_6) +
    stat_poly_line(formula = formula_6)+
    theme_classic()
ggsave(file.path(out_dir, "effect_size_vs_library_size.png"), p6, height = 3, width = 3)

# Make plot: library size vs proportion of intron detected
formula_7 <- y ~ x + poly(x, 2)
p7 <- ggplot(cell_type_summary_2, aes(x = mean_libsize, y = mean_prop_detected)) +
    geom_point() +
    theme_classic() +
    xlab("Library size") +
    ylab("Proportion of detected introns") +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_7) +
    stat_poly_line(formula = formula_7)
ggsave(file.path(out_dir, "prop_detected_introns_vs_library_size.png"), p7, height = 3, width = 3)

# Is library size and sample size independent contributors?
fit1 <- lm(formula = mean_slope ~ n_donors + poly(n_donors, 2), cell_type_summary_2)
cell_type_summary_3 <- cell_type_summary_2 %>%
    mutate(
        residual_slope = residuals(fit1),
        pred_slope = predict(fit1)
    )

# Make plot: residual effect size vs library size
formula_8 <- y ~ x
p8 <- ggplot(cell_type_summary_3, aes(x = mean_libsize, y = residual_slope)) +
    geom_point() +
    xlab("Library size") +
    ylab("Mean absolute effect size") +
    stat_poly_eq(use_label(c("adj.R2")), formula = formula_8, label.x = "right", label.y = "bottom") +
    stat_poly_line(formula = formula_8) +
    theme_classic()
ggsave(file.path(out_dir, "residual_effect_size_vs_library_size.png"), p8, height = 3, width = 3)


fit <- lm(formula = sgene ~ n_donors + mean_libsize, data = cell_type_summary_1)
anova(fit)
#              Df  Sum Sq Mean Sq F value    Pr(>F)
# n_donors      1 2296537 2296537  256.72 1.716e-12 ***
# mean_libsize  1 2366134 2366134  264.50 1.316e-12 ***
# Residuals    19  169970    8946
