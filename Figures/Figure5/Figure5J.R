library(dplyr)
library(ggplot2)
library(data.table)
library(readr)
library(ggnewscale)
library(cowplot)

out_pdf_path = "/ebs1/users/liufei/project/tchp20230625/data/processed/snp_2_finemapping_3celltype_20230712/plot_combine_pip_logP_eqtl_em_naive_v2.pdf"

merge_path <- "/ebs1/users/jingzhixuan/AIDA_sQTL/HNRNPLL-PTPRC-coloc/locus/"
ld_path <- "/ebs1/users/jingzhixuan/AIDA_sQTL/HNRNPLL-PTPRC-coloc/locuscompare/"
snp_naive <- "rs6751481"
snp_em <- "rs74258942"

naive_eqtl <- fread(paste(merge_path,"CD4+_T_naive.txt",sep=""))
em_eqtl <- fread(paste(merge_path,"CD4+_T_em.txt",sep=""))

ld_naive <- fread(paste(ld_path,"CD4+_T_naive/ld.tsv",sep=""))
ld_em <- fread(paste(ld_path,"CD4+_T_em/ld.tsv",sep=""))
ld_naive <- rbind(ld_naive, data.frame(SNP_A = snp_naive, SNP_B = snp_naive, R2 = 1))
ld_em <- rbind(ld_em, data.frame(SNP_A = snp_em, SNP_B = snp_em, R2 = 1))

naive_eqtl <- naive_eqtl %>% select(rsid, pos, epval, logpe)
naive_eqtl <- naive_eqtl[(!duplicated(naive_eqtl[, c("pos")])), ]
naive_eqtl <- left_join(naive_eqtl, ld_naive, by = c("rsid" = "SNP_B"))
naive_eqtl$celltype <- "CD4+_T_naive"

em_eqtl <- em_eqtl %>% select(rsid, pos, epval, logpe)
em_eqtl <- em_eqtl[(!duplicated(em_eqtl[, c("pos")])), ]
em_eqtl <- left_join(em_eqtl, ld_em, by = c("rsid" = "SNP_B"))
em_eqtl$celltype <- "CD4+_T_em"

naive_eqtl$R2 <- ifelse(is.na(naive_eqtl$R2), 0.00001, naive_eqtl$R2)
em_eqtl$R2 <- ifelse(is.na(em_eqtl$R2), 0.00001, em_eqtl$R2)
shape_naive <- ifelse(naive_eqtl$rsid == snp_naive, 15, 17)
shape_em <- ifelse(em_eqtl$rsid == snp_em, 15, 16)
size_naive <- ifelse(naive_eqtl$rsid == snp_naive, 3, 2)
size_em <- ifelse(em_eqtl$rsid == snp_em, 3, 2)
title <- "cis-eQTL"

eqtl <- ggplot(mapping = aes(x = pos, y = logpe)) +
    geom_point(data = naive_eqtl, aes(color = R2), shape = shape_naive, size = size_naive) +
    scale_color_gradientn(colors = c('#FFCCCC', '#FF0000')) +
    new_scale_color() +
    geom_point(data = em_eqtl, aes(color = R2, x = pos, y = logpe*10/11), shape = shape_em, size = size_em) +
    theme_classic()+
    scale_color_gradientn(colors = c('#CCCCFF', '#0000FF')) +
    xlim(38600000, 38800000) +
    # scale_y_continuous(name = "logp-CD4+T naive", expand = c(0,0),limits = c(0,8), sec.axis = sec_axis(trans=~.*11/10, name="logp-CD4+T em")) +
    scale_y_continuous(expand = c(0,0),limits = c(0,8), sec.axis = sec_axis(trans=~.*11/10, name=expression(paste("em CD4+T cis-eQTL-", log[10],"(P)",sep="")))) +
    
    xlab('chr 2 (bp)') +
    # ylab(bquote(.(title)~-log[10]*'(P)')) +
    # ylab(expression(paste("Naive CD4+T cis-eQTL-log10(P)"))) +
    ylab(expression(paste("Naive CD4+T cis-eQTL-", log[10],"(P)",sep="")))+
    theme(
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),

        axis.line.y.right = element_line(color = "#0000FF"),
        axis.line.y.left = element_line(color = "#FF0000"),
        # axis.ticks.y.right = element_line(color = "#0000FF"),
        # axis.text.y.right = element_text(color = "#0000FF"),
        # axis.title.y.right = element_text(color = "#0000FF"),
        # legend.position="none"
    )
    




naive_pip = read.csv("/ebs1/users/liufei/project/tchp20230625/data/processed/snp_2_finemapping_3celltype_20230712/CD4+_T_naive.eQTL.full.rsid.susiepip.tsv",sep='\t')
em_pip = read.csv("/ebs1/users/liufei/project/tchp20230625/data/processed/snp_2_finemapping_3celltype_20230712/CD4+_T_em.eQTL.full.rsid.susiepip.tsv",sep='\t')
em_pip$result.pip<-em_pip$result.pip*3


naive_susiepip <- left_join(naive_pip, ld_naive, by = c("rsid" = "SNP_B"))
naive_susiepip$celltype <- "CD4+_T_naive"

em_susiepip <- left_join(em_pip, ld_em, by = c("rsid" = "SNP_B"))
em_susiepip$celltype <- "CD4+_T_em"

naive_susiepip$R2 <- ifelse(is.na(naive_susiepip$R2), 0.00001, naive_susiepip$R2)

em_susiepip$R2 <- ifelse(is.na(em_susiepip$R2), 0.00001, em_susiepip$R2)

shape_naive <- ifelse(naive_susiepip$rsid == snp_naive, 15, 17)

shape_em <- ifelse(em_susiepip$rsid == snp_em, 15, 16)

size_naive <- ifelse(naive_susiepip$rsid == snp_naive, 3, 2)

size_em <- ifelse(em_susiepip$rsid == snp_em, 3, 2)



pip <- ggplot() +
    geom_point(data=naive_susiepip, aes(x = pos, y = result.pip, color = R2), shape = shape_naive, size = size_naive) +
    scale_color_gradientn(colors = c('#FFCCCC', '#FF0000')) +
    new_scale_color() +
    
    geom_point(data=em_susiepip, aes(x = pos, y = result.pip, color = R2), shape = shape_em, size = size_em) +
    scale_color_gradientn(colors = c('#CCCCFF', '#0000FF')) +
    xlab('chr 2 (bp)') +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . *1/3, name = "em CD4+T cis-eQTL PIP",breaks = c(0, 0.075, 0.15, 0.225))) +
    theme_classic() +

    theme(
        # axis.text.x = element_blank(), # no label
        # axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),

    #     #legend.position = c(0.95, 0.95),
    #     #legend.justification = c("right", "top"),
    #     #legend.background = element_rect(
    #     #  fill = alpha("grey90", 0.8),
    #     #  color = "black",
    #     #  linewidth = 0.5, linetype = "solid"
    #     #),
        axis.line.y.right = element_line(color = "#0000FF"),
        axis.line.y.left = element_line(color = "#FF0000"),
        # axis.ticks.y.right = element_line(color = "#0000FF"),
        # axis.text.y.right = element_text(color = "#0000FF"),
        # axis.title.y.right = element_text(color = "#0000FF"),
        # legend.position="none"
    )+
    ylab("Naive CD4+T cis-eQTL PIP") 


    final <- plot_grid(eqtl + theme(legend.position="none"),
                      pip + theme(legend.position="none"),
                      ncol = 1,
                      align = "v") +
            geom_vline(xintercept = c(38670668, 38706826), linetype = "dotted")
    # create some space to the left of the legend
    legend <- get_legend(pip)

    final_legend <- plot_grid(final, legend, rel_widths = c(3, .4))

    ggsave(out_pdf_path, final_legend, width = 6, height = 4)





