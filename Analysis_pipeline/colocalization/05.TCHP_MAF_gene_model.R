###This script was used to generate TCHP's allele frequency both in 1000 Genomes and AIDA dataset. Gene model of how allele mutation influences TCHP
###splicing was also shown here.

library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)

maf_file <- "/ebs1/users/boxiangliu/projects/aida_sqtl/scripts/exp_validation/AF_original_data.tsv"
fig_file <- "/ebs1/users/boxiangliu/projects/aida_sqtl/scripts/exp_validation/AF_plot.png"

maf <- fread(maf_file)
maf <- maf %>% mutate(MAF = ifelse(MAF == 0, -0.0001, MAF))
maf <- maf %>%
    mutate(EAS = ifelse(population %in% c("Chinese", "Japanese", "Korean", "EAS"), "EAS", "Others"))

order1 <- maf %>%
    filter(data == "1000 Genomes") %>%
    arrange(desc(MAF)) %>%
    select(population) %>%
    unlist() %>%
    unname()

order2 <- maf %>%
    filter(data == "AIDA") %>%
    arrange(desc(MAF)) %>%
    select(population) %>%
    unlist() %>%
    unname()

p1 <- ggplot(maf[data == "1000 Genomes"], aes(x = population, y = MAF, fill = EAS)) +
    geom_bar(stat = "identity", color = "black") +
    scale_x_discrete(limits = order1) +
    scale_fill_manual(values = c("darkolivegreen3", "black")) +
    theme_classic() +
    xlab(NULL) +
    theme(legend.position = c(0.85, 0.55), legend.background = element_rect(fill = "transparent"), legend.title = element_blank()) +
    annotate(geom = "text", label = "1000 Genomes", x = 4.2, y = 0.08)

p2 <- ggplot(maf[data == "AIDA"], aes(x = population, y = MAF, fill = EAS)) +
    geom_bar(stat = "identity", color = "black") +
    theme_classic() +
    scale_x_discrete(limits = order2) +
    scale_fill_manual(values = c("darkolivegreen3", "black")) +
    xlab(NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    theme(legend.position = "none") +
    annotate(geom = "text", label = "AIDA", x = 5, y = 0.15)

p3 <- p2 / p1
ggsave(fig_file, p3, height = 3, width = 3)


###Gene model of how allele mutation influences TCHP splicing
library(ggplot2)
library(gggenes)
library(data.table)
library(dplyr)
library(ggsci)

tchp <- fread("/ebs1/users/boxiangliu/projects/aida_sqtl/scripts/exp_validation/TCHP_model.tsv")
features <- fread("/ebs1/users/boxiangliu/projects/aida_sqtl/scripts/exp_validation/TCHP_features.tsv")
plot_file <- "/ebs1/users/boxiangliu/projects/aida_sqtl/scripts/exp_validation/TCHP_model.png"

tchp <- tchp %>% mutate(exon4 = start == 109904737)

p <- ggplot(tchp, aes(xmin = start, xmax = end, y = transcript, fill = exon4, color = exon4)) +
    geom_gene_arrow(arrowhead_width = grid::unit(0, "mm"), arrowhead_height = grid::unit(3, "mm")) +
    scale_fill_manual(values = c("black", "red")) +
    scale_color_manual(values = c("black", "red")) +
    theme_genes() +
    geom_feature(
        data = features,
        aes(x = position, y = transcript, forward = forward)
    ) +
    xlab("chr12") +
    ylab(NULL) +
    theme(legend.position = "none")
ggsave(plot_file, p, width = 5, height = 1.5)
