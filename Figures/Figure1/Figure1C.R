library(ggplot2)
library(ggExtra)
library(dplyr)

# read in 
dat1 <- read.delim('data/mergeAll001.all_celltypes.R1.exon.cov.bed', header = F)
dat2 <- read.delim('data/mergeAll001.all_celltypes.R2.exon.cov.bed', header = F)

## filter overlapped genes
overlap_genelist1 <- unlist(sapply(unique(dat1$V5[grep(',', dat1$V5)]), function(x)strsplit(x, ',')))
overlap_genelist2 <- unlist(sapply(unique(dat2$V5[grep(',', dat2$V5)]), function(x)strsplit(x, ',')))
length(overlap_genelist1)
length(overlap_genelist2)
overlap_genelist <- unique(overlap_genelist1)
length(overlap_genelist)

dat1.nonoverlap <- subset(dat1, ! V5 %in% c(overlap_genelist, unique(dat1$V5[grep(',', dat1$V5)])))
dat2.nonoverlap <- subset(dat2, ! V5 %in% c(overlap_genelist, unique(dat2$V5[grep(',', dat2$V5)])))

# calculate perc covered 
dat1.nonoverlap.group <- dat1.nonoverlap %>%
  group_by(V5) %>%
  summarise(all_bases = sum(V9), covered_bases = sum(V8))
dat1.nonoverlap.group$perc_covered_bases <- dat1.nonoverlap.group$covered_bases/dat1.nonoverlap.group$all_bases
dat2.nonoverlap.group <- dat2.nonoverlap %>%
  group_by(V5) %>%
  summarise(all_bases = sum(V9), covered_bases = sum(V8))
dat2.nonoverlap.group$perc_covered_bases <- dat2.nonoverlap.group$covered_bases/dat2.nonoverlap.group$all_bases


# filter expressed genes
expr_genelist1 <- dat1 %>%
  group_by(V5) %>%
  summarise(sumReads = sum(V7))
expr_genelist2 <- dat2 %>%
  group_by(V5) %>%
  summarise(sumReads = sum(V7))
#dat.nonoverlap.group.expr <- subset(dat.nonoverlap.group, V5 %in% subset(expr_genelist, sumReads > 5)$V5 & covered_bases >= 112)
dat1.nonoverlap.group <- merge(dat1.nonoverlap.group, expr_genelist1)
dat2.nonoverlap.group <- merge(dat2.nonoverlap.group, expr_genelist2)
dat1.nonoverlap.group.expr <- subset(dat1.nonoverlap.group, V5 %in% subset(expr_genelist1, sumReads > 5)$V5 & V5 %in% subset(expr_genelist2, sumReads > 5)$V5)
dat2.nonoverlap.group.expr <- subset(dat2.nonoverlap.group, V5 %in% subset(expr_genelist1, sumReads > 5)$V5 & V5 %in% subset(expr_genelist2, sumReads > 5)$V5)

dat1.nonoverlap.group.expr$read <- 'R1'
dat2.nonoverlap.group.expr$read <- 'R2'
dat.nonoverlap.group.expr <- rbind(dat1.nonoverlap.group.expr, dat2.nonoverlap.group.expr)
median(dat.nonoverlap.group.expr$perc_covered_bases)
mean(dat.nonoverlap.group.expr$perc_covered_bases)







# merged R1 R2
dat_merge <- read.delim('data/mergeAll001.all_celltypes.merge.R1_R2rev.exon.cov.bed', header = F)
## filter overlapped genes
dat_merge.nonoverlap <- subset(dat_merge, ! V5 %in% c(overlap_genelist, unique(dat_merge$V5[grep(',', dat_merge$V5)])))

# calculate perc covered 
dat_merge.nonoverlap.group <- dat_merge.nonoverlap %>%
  group_by(V5) %>%
  summarise(all_bases = sum(V9), covered_bases = sum(V8))
dat_merge.nonoverlap.group$perc_covered_bases <- dat_merge.nonoverlap.group$covered_bases/dat_merge.nonoverlap.group$all_bases


# filter expressed genes
expr_genelist <- dat_merge %>%
  group_by(V5) %>%
  summarise(sumReads = sum(V7))

dat_merge.nonoverlap.group <- merge(dat_merge.nonoverlap.group, expr_genelist)
dat_merge.nonoverlap.group.expr <- subset(dat_merge.nonoverlap.group, V5 %in% subset(expr_genelist1, sumReads > 5)$V5 & V5 %in% subset(expr_genelist2, sumReads > 5)$V5)
median(dat_merge.nonoverlap.group.expr$perc_covered_bases)
# 0.9452289
mean(dat_merge.nonoverlap.group.expr$perc_covered_bases)
# 0.8081093


panel1 <- ggplot(dat_merge.nonoverlap.group.expr) + geom_boxplot(outlier.size=.1,aes(x=bin,y=perc_covered_bases)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Exonic length (bp)') + ylab('') 
panel2 <- ggplot(dat_merge.nonoverlap.group.expr) + geom_point(size=.1,alpha=0.4,aes(x=sumReads,y=perc_covered_bases)) + theme_bw() + 
  xlab('# of reads') + ylab('') + scale_x_log10() #+ geom_smooth(method = 'lm', se=TRUE, formula = y ~ x, orientation = 'y', aes(x=sumReads,y=perc_covered_bases)) 


#write.table(dat_merge.nonoverlap.group.expr$V5, 'genelist_expressed.txt', quote = F, row.names = F, col.names = F)
library(grid)
library(gridExtra)
# 2 panel 
dat_merge.nonoverlap.group.expr = dat_merge.nonoverlap.group.expr %>% mutate(
  bin_numReads = case_when(
    sumReads <= 10 ~ "1 - 10",
    sumReads <= 100 ~ "10 - 100",
    sumReads <= 1000 ~ "100 - 1000",
    sumReads <= 10000 ~ "1000 - 10000",
    sumReads <= 100000 ~ "10000 - 100000",
    TRUE ~ "100000+" 
  )
)
table(dat_merge.nonoverlap.group.expr$bin_numReads)
df_bin_numReads <- as.data.frame(table(dat_merge.nonoverlap.group.expr$bin_numReads))
df_bin_numReads <- subset(df_bin_numReads, Freq != 0)
dat_merge.nonoverlap.group.expr$bin_numReads <- factor(dat_merge.nonoverlap.group.expr$bin_numReads, 
                                                levels = c("1 - 10",
                                                         "10 - 100",
                                                         "100 - 1000",
                                                         "1000 - 10000",
                                                         "10000 - 100000",
                                                         "100000+"))


panel2a <- ggplot(dat_merge.nonoverlap.group.expr) + geom_boxplot(outlier.size=.1, fatten=1, aes(x=bin_numReads,y=perc_covered_bases, fill=bin_numReads)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab('') + ylab('')  + scale_fill_discrete(name = "# of reads per gene")
panel2aa <- ggplot(dat_merge.nonoverlap.group.expr) + geom_boxplot(outlier.shape = NA, fatten=1, aes(x=bin_numReads,y=perc_covered_bases, fill=bin_numReads)) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab('') + ylab('Percent of exonic bases covered')  + 
  scale_fill_manual(name = "# of reads per gene", values = colorRampPalette(c("white", "#686868"))(5)) +
  geom_text(inherit.aes = F, aes(x = Var1, y = .1, label = Freq), df_bin_numReads)
panel3 <- ggplot(dat.nonoverlap.group.expr, aes(x = '', y = perc_covered_bases, fill = '#7A7F3D')) + theme_classic() +
  geom_violin(alpha = .4) + xlab('') +
  geom_boxplot(width = .1, alpha = .6, outlier.shape = NA) +
  scale_fill_manual(values = c('#7A7F3D')) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ylab('') + geom_hline(yintercept = mean(dat.nonoverlap.group.expr$perc_covered_bases), linetype = 2, color = 'red')

quantile(dat.nonoverlap.group.expr$perc_covered_bases)
#0%        25%        50%        75%       100% 
#0.00543418 0.51805786 0.85346447 0.97267230 1.00000000 

library(cowplot)
panel23a <- plot_grid(panel2aa, panel3, ncol = 2, rel_widths = c(3/5,2/5))
ggsave('Fig.1c.pdf', panel23a, width = 8, height = 3)


