library(reshape2)
celltypes <- read.table('data/cell_type_use.QTL.txt')$V1
dat <- read.delim('data/stats_pi1.cell_type.tsv')
dat <- subset(dat, QTL_celltype %in% celltypes & pi1_celltype %in% celltypes)
dat.cis <- read.delim('data/pi1_cell_type_cis_sQTL', check.names = F)

dat$id <- paste(dat$QTL_celltype, dat$pi1_celltype, sep = '__')

dat.cis$ct_r <- rownames(dat.cis)
dat.cis.long <- melt(dat.cis, id.vars = 'ct_r', variable.name = 'ct_c', value.name = 'pi1_cis')
dat.cis.long <- subset(dat.cis.long, ct_r != ct_c)
dat.cis.long$id <- paste(dat.cis.long$ct_r, dat.cis.long$ct_c, sep = '__')

dat.merge <- merge(dat, dat.cis.long)
dat.merge.filter <- na.omit(dat.merge)
dat.merge.filter.melt <- melt(dat.merge.filter, measure.vars = c('pi1', 'pi1_cis'), id.vars = 'id', value.name = 'pi1', variable.name = 'QTL')
library(ggplot2)
library(ggpubr)
library(dplyr)
dat.merge.filter.melt <- dat.merge.filter.melt %>%
  mutate(QTL = case_when(QTL == 'pi1' ~ 'trans', QTL == 'pi1_cis' ~ 'cis'))
ggpaired(dat.merge.filter.melt, x = "QTL", y = "pi1",
         color = "QTL", line.color = "gray", line.size = 0.1,
         palette = "npg",
         xlab = 'QTL', ylab = 'pi1') + 
  stat_compare_means(paired = TRUE, method = 'wilcox.test') 
ggsave('Fig.5c.pdf', width = 3.5, height = 3.5)
