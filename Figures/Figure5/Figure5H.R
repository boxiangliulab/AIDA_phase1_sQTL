library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(lme4)
celltypes <- read.table('data/cell_type_use.txt')$V1
ancestry_PC <- read.delim('data/table.ancestry_PC.tsv')
meta <- read.csv('data/age_sex.csv')
# read files
geno <- read.delim('data/geno.txt')
pheno <- read.delim('data/num_cells_per_donor.tsv')
pheno.mat <- dcast(pheno, DCP_ID ~ Annotation_Level2, value.var = 'n')
pheno.mat$naive_over_mature <- pheno.mat$`CD4+_T_naive`/(pheno.mat$`CD4+_T_cm`+pheno.mat$`CD4+_T_em`+pheno.mat$`CD4+_T_cyt`)
pheno.mat$log10_naive_over_mature <- log10(pheno.mat$naive_over_mature)
# preprocessing
i <- 6  ##change according to your need
#i <- 7  ##change according to your need

geno.fmt <- as.data.frame(t(geno[i, grep('H[0-9]', colnames(geno))]))
colnames(geno.fmt) <- 'genotype'
geno.fmt <- geno.fmt %>%
  separate(genotype, unlist(strsplit(geno$FORMAT[i], split = ':')), sep = ':')
table(geno.fmt$GT)
geno.num <- geno.fmt


# 
ref<-geno$REF[i]
alt<-geno$ALT[i]

# replace 0:0 0:1 1:0 1:1
for(k in 1:nrow(geno.fmt)){
  if(startsWith(geno.fmt[k,1],'0|0'))geno.fmt[k,1]<-paste(ref,ref,sep="")
  if(startsWith(geno.fmt[k,1],'0|1'))geno.fmt[k,1]<-paste(ref,alt,sep="")
  if(startsWith(geno.fmt[k,1],'1|0'))geno.fmt[k,1]<-paste(ref,alt,sep="")
  if(startsWith(geno.fmt[k,1],'1|1'))geno.fmt[k,1]<-paste(alt,alt,sep="")
}
for(k in 1:nrow(geno.num)){
  if(startsWith(geno.num[k,1],'0|0'))geno.num[k,1]<-0
  if(startsWith(geno.num[k,1],'0|1'))geno.num[k,1]<-1
  if(startsWith(geno.num[k,1],'1|0'))geno.num[k,1]<-1
  if(startsWith(geno.num[k,1],'1|1'))geno.num[k,1]<-2
}
pheno$GT <- factor(geno.fmt[pheno$DCP_ID, 'GT'], 
                   levels = c(paste(ref,ref,sep=""), paste(ref,alt,sep=""), paste(alt,alt,sep="")))
pheno$alt_freq <- geno.num[pheno$DCP_ID, 'GT']


pheno.mat$GT <- factor(geno.fmt[pheno.mat$DCP_ID, 'GT'],
                       levels = c(paste(ref,ref,sep=""), paste(ref,alt,sep=""), paste(alt,alt,sep="")))
pheno.mat$alt_freq <- geno.num[pheno.mat$DCP_ID, 'GT']

## linear model
pheno.mat.lm <- na.omit(pheno.mat[, c('DCP_ID', 'alt_freq', 'naive_over_mature', 'log10_naive_over_mature')])
pheno.mat.lm$alt_freq <- as.numeric(pheno.mat.lm$alt_freq)
# genotype PC
pheno.mat.lm <- merge(pheno.mat.lm, ancestry_PC, by.x = 'DCP_ID', by.y = 'ident')
# sequence center
pheno.mat.lm$seq_center <- substr(pheno.mat.lm$DCP_ID, 1, 6)
pheno.mat.lm <- pheno.mat.lm %>% mutate(dummy=1) %>% pivot_wider(names_from = seq_center, values_from = dummy, values_fill = 0)
# age, sex
pheno.mat.lm <- merge(pheno.mat.lm, meta, by = 'DCP_ID')


fm <- lm(log10_naive_over_mature ~ alt_freq + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + age + sex + JP_RIK + KR_SGI ,data=pheno.mat.lm)
coefs <- data.frame(coef(summary(fm)))
coefs$p.value <- 2 * pnorm(abs(coefs$t.value), lower.tail = F)



pheno.mat.lm$alt_freq <- gsub('0', paste0(ref,ref), pheno.mat.lm$alt_freq)
pheno.mat.lm$alt_freq <- gsub('1', paste0(ref,alt), pheno.mat.lm$alt_freq)
pheno.mat.lm$alt_freq <- gsub('2', paste0(alt,alt), pheno.mat.lm$alt_freq)








library(dplyr)
library(reshape2)
total.all.sig.endpoints.group_mean <- pheno.mat.lm %>%
  group_by(alt_freq) %>%
  summarise(mean = mean(log10_naive_over_mature)) %>%
  dcast(0~alt_freq, value.var = 'mean')

total.all.sig.endpoints.group_mean$x1 <- 'TT'
total.all.sig.endpoints.group_mean$x2 <- 'CC'

pheno.mat.lm$alt_freq <- factor(pheno.mat.lm$alt_freq, levels = c(paste0(ref,ref), paste0(ref,alt), paste0(alt,alt)))
ggplot(pheno.mat.lm, aes(x = alt_freq, y = log10_naive_over_mature)) +
  geom_violin(color = '#807dba', alpha = .7) +
  geom_boxplot(color = '#807dba', width=0.05) + theme_classic() +
  ylab('Log10(Naive CD4+ T/(cm CD4+ T + em CD4+ T)') + xlab('Genotype)') +
  ggtitle(paste0(geno$ID[i], 
                 '\nP-value: ', signif(coefs$p.value[2], 2),
                 '\nBeta: ', signif(coefs$Estimate[2], 2))) +
  geom_segment(aes(x = x1, y = TT, xend = x2, yend = CC), color = 'red', data = total.all.sig.endpoints.group_mean) 


ggsave('Fig.5h.pdf', width = 3, height = 3.5)

table(pheno.mat.lm$alt_freq)
