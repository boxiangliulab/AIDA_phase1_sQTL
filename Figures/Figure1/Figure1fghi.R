#' this is to plot the newest version of figure1 F, combining splicing result from AIDA (5') with onek1k (3')
#' date: 2023.06.21
aida<-read.table("AIDA_fist_phase_metadata_add_annotation_add_spliz_result_0618_subset.tsv",sep="\t",header = T)
aida<-dplyr::filter(aida,celltype!="CD4+_T")

## change name 
label_dict<-list("CD14\\+_Monocyte"="CD14+ Monocyte",
                 "CD16\\+_Monocyte"="CD16+ Monocyte",
                 "cDC2"="cDC2",
                 "CD16\\+_NK"="CD16+ NK",
                 "CD56\\+_NK"="CD56+ NK",
                 "CD4\\+_T_cm"="cm CD4+ T",
                 "CD4\\+_T_cyt"="cyt CD4+ T",
                 "CD4\\+_T_em"="em CD4+ T",
                 "CD4\\+_T_naive"="Naive CD4+ T",
                 "Treg"="Treg",
                 "CD8\\+_T_naive"="Naive CD8+ T",
                 "CD8\\+_T_GZMB\\+"="GZMBhi CD8+ T",
                 "CD8\\+_T_GZMK\\+"="GZMKhi CD8+ T",
                 "MAIT"="MAIT",
                 "gdT_GZMK\\-"="GZMBhi gdT",
                 "gdT_GZMK\\+"="GZMKhi gdT",
                 "naive_B"="Naive B",
                 "IGHMlo_memory_B"="IGHMlo memory B",
                 "IGHMhi_memory_B"="IGHMhi memory B",
                 "pDC"="pDC",
                 "atypical_B"="Atypical B",
                 "CD4\\+_T" = "CD4+ T") # 22 cell types

#' replace old names from the plot file
#  plot_data is your dataframe
for (i in seq_along(label_dict)) {
  aida$celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],aida$celltype)
}

colors2<-c(
  "CD14+ Monocyte"="#fee090", "CD16+ Monocyte"="#fdae61", "CD16+ NK"="#bf812d",
  "CD4+ T"="#4a1486", "cm CD4+ T"="#9e9ac8", "cyt CD4+ T"="#6a51a3", "em CD4+ T"="#807dba",
  "Naive CD4+ T"="#bcbddc", "CD56+ NK"="#8c510a", "GZMBhi CD8+ T"="#3690c0", "GZMKhi CD8+ T"="#74a9cf",
  "Naive CD8+ T"="#a6bddb", "IGHMhi memory B"="#a1d99b", "IGHMlo memory B"="#74c476",
  "MAIT"="#0570b0", "Treg"="#CCCC4D", "Atypical B"="#238b45", "cDC2"="#f46d43",
  "GZMKhi gdT"="#f768a1", "GZMBhi gdT"="#dd3497", "Naive B"="#c7e9c0", "pDC"="#d73027")



#order number of spliced genes
aida$Annotation_Level1[grep("DC",aida$Annotation_Level1)]<-"DCs"
aida$Annotation_Level1<-factor(aida$Annotation_Level1,levels = c("DCs","Monocyte","B","CD4+_T","CD8+_T","gdT","NK"))
summary <- aida %>%
  group_by(Annotation_Level1,celltype) %>%
  summarise(median_n = median(Num_genes))

summary<-summary %>%  group_by(Annotation_Level1) %>% arrange(desc(median_n),.by_group = TRUE)

aida$celltype <- factor(aida$celltype, levels = summary$celltype)

# caclulate average nFeature

nfeature<-aida %>%
  group_by(celltype) %>%
  summarise(mean_nFeature = mean(nFeature_RNA))

nfeature$celltype<-factor(nfeature$celltype,levels = summary$celltype)
#nfeature$mean_nFeature<-nfeature$mean_nFeature*(0.75) scale

# caclulate average ratio

aida$spliz_ratio<-aida$Num_genes/aida$nFeature_RNA
spliz_ratio<-aida %>% 
  group_by(celltype) %>%
  summarise(mean_spliz_ratio = mean(spliz_ratio))
spliz_ratio$celltype<-factor(spliz_ratio$celltype,levels = summary$celltype)
spliz_ratio$mean_spliz_ratio<-spliz_ratio$mean_spliz_ratio*5000
pdf("Fig1F.spliz_with_onek1k.v6.pdf",width = 5 ,height = 3)
p1<-ggplot() +
  geom_boxplot(data = aida, aes(x = celltype, y = Num_genes,fill=celltype),size=0.1,outlier.size=0.5) +
  theme_classic() +
  scale_fill_manual(values=colors2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_point(data = nfeature, aes(x = celltype, y = mean_nFeature), size = 2, color = "black",fill="red",shape=23,stroke=0.2) +
  scale_x_discrete(name = "celltype")+
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))+
  geom_hline(yintercept = 267,linetype="dashed",color="#4DBBD5FF")+
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position="none"
  )+
  ylab("Number of AS genes")

dev.off()

# plot distribution
# 3' 12 individual distribution with aida

count_12<-read.table("12indiv_after_filter_merge_expression_splice.tsv",header = T)
count<-data.frame(Num_genes=count_12$Num_genes)
count$type="Onek1k"
#extract spliced gene num from AIDA
aida_count<-data.frame(Num_genes=aida$Num_genes)
aida_count$type="AIDA"
count_merge<-rbind(aida_count,count)

## plot aida and onek1k together

set.seed(123)
aida.extr <- aida[sample(nrow(aida),nrow(count_12)),c(3,4)] %>% mutate(type="AIDA") # sample the number of cells that equal to onek1k
#aida.extr<-aida[,c(3,4)] %>% mutate(type="AIDA")
onek1k.extr<-count_12[,c(4,7)] %>% mutate(type="Onek1k")

equition.aida<-paste("y =", round(coef(summary(lm(Num_genes ~ nFeature_RNA , data = aida.extr)))[2, "Estimate"], 2), "x", round(coef(summary(lm(Num_genes ~ nFeature_RNA , data = aida.extr)))[1, "Estimate"], 2))
equition.onek1k<-paste("y =", round(coef(summary(lm(Num_genes ~ nFeature_RNA , data = onek1k.extr)))[2, "Estimate"], 2), "x +", round(coef(summary(lm(Num_genes ~ nFeature_RNA, data = onek1k.extr)))[1, "Estimate"], 2))
merge<-rbind(aida.extr,onek1k.extr)
p2<-ggplot(merge, aes(x = nFeature_RNA, y = Num_genes, color = type)) +
  geom_point(size = 0.25) +
  scale_color_manual(values = c("AIDA" = "#E64B35FF", "Onek1k" = "#4DBBD5FF")) +
  geom_smooth(method = "lm", se = FALSE, aes(group = type), color = "black", linewidth = 0.3) +
  theme_classic() +
  ylab("Detected spliced gene number") +
  xlab("Detected expressed gene number") +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(color = "") +
  theme(
    legend.justification = c(0.05, 0.95),
    legend.position = c(0.05, 0.95),
    #legend.box.background = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

## calculate confidential interval when p-value=0.05

mean(aida$spliz_ratio)
#[1] 0.5775968

aida.ci<-t.test(aida$spliz_ratio,conf.level = 0.95)


onek1k.ci<-t.test(count_12$splice_ratio,conf.level = 0.95)
mean(count_12$splice_ratio)

######## ols
##### aida subset

x<-aida.extr$nFeature_RNA
y=aida.extr$Num_genes
fit_aida<-lm(y~x)
aida.ci.ols<-confint(fit_aida, level = 0.95)

# ### test linear relationship aida

# 12 individual
x<-aida$nFeature_RNA
y<-aida$Num_genes
fit<-lm(y~x)
summary(fit)

##plot figure h
data<-read.csv("gene_donor_num.csv")
data<-dplyr::filter(data,Celltype!="CD4+_T")
dat.filter<-data
for (i in seq_along(label_dict)) {
  dat.filter$Celltype <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],dat.filter$Celltype)
}
summary <- dat.filter %>% group_by(Celltype) %>% summarise(median_n = median(num))
dat.filter$Celltype <- factor(dat.filter$Celltype, levels = summary$Celltype[order(summary$median_n, decreasing = T)])

pdf("leafcutter_plot_gene_per_cell_v6.pdf",width = 4,height = 3)

p3<-ggplot(dat.filter,aes(x=Celltype,y=num,fill=Celltype)) +
  geom_boxplot(size=0.1,outlier.size=0.5)+ 
  scale_fill_manual(values=colors2)+theme_classic() +xlab("Cell type") + 
  ylab("# of detected genes")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 10),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))+
  ylab("Number of AS genes")

dev.off()

##plot figure i
library(readxl)
library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggpmisc)

cell_num<-read.table("/Users/tongyihan/project/AIDA/data/figure/figure1/data/22_cell_type_individual_over_10_cells_no_flag.csv",sep=",",header=TRUE)
gene_num<-read.csv("/Users/tongyihan/project/AIDA/data/figure/figure1/gene_donor_num.csv")

for(i in 1:8190){
  cell_num[i,2]<-paste(cell_num[i,1],cell_num[i,2],sep="")
  gene_num[i,2]<-paste(cell_num[i,1],gene_num[i,2],sep="")
}
colnames(cell_num)[2]<-"Celltype"

newdata<-merge(cell_num,gene_num,by="Celltype")

## fit nls model
fit_nls <- nls(num ~ SSlogis(log10(number), Asym, xmid, scal), newdata)


fit_summary <- summary(fit_nls)
formula <- formula(fit_nls)
parameters <- coef(fit_nls)

print(fit_summary)
print(formula)
print(parameters)
a <- parameters["Asym"] # a numeric parameter representing the asymptote.
b <- parameters["xmid"] # a numeric parameter representing the x value at the inflection point of the curve. The value of SSlogis will be Asym/2 at xmid.
predictions <- predict(fit_nls, newdata)
#formula_text <- paste0("y =", round(coef(fit_nls)['Asym'], 2), " / (1 + exp(-", round(coef(fit_nls)['xmid'], 2), " * (log10(x) - ", round(coef(fit_nls)['scal'], 2), ")))")
#formula_text <- paste0("y =", round(coef(fit_nls)['Asym'], 2), " / (1 + exp(", round(coef(fit_nls)['xmid'], 2), " * (log10(x) - ", round(coef(fit_nls)['scal'], 2), ")))")

formula_text<-paste0("y =", round(coef(fit_nls)['Asym'], 2), " / (1 + exp((1.21 - log10(x)) / ", round(coef(fit_nls)['scal'], 2), "))")
formula_text<-paste0("y =", round(coef(fit_nls)['Asym'], 2), " / (1 + exp((1.21 - x) / ", round(coef(fit_nls)['scal'], 2), "))")

## calculate R2

predicted <- predict(fit_nls, newdata)

TSS <- sum((newdata$num - mean(newdata$num))^2)

RSS <- sum((newdata$num - predicted)^2)

R_squared <- 1 - (RSS / TSS)

print(R_squared)
#[1] 0.9265301

#my.formula <- y ~ log10(x)

p4<-ggplot(newdata, aes(x = number, y = num)) +
  geom_point(size = 0.3,color="#E64B35FF") +
  geom_line(data = data.frame(number = newdata$number, num = predictions), color = "black") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(
    axis.text.x = element_text(size = 14),
    #axis.title.x = element_text(size = 16,hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.justification = c(0.05, 0.95),
    legend.position = c(0.05, 0.95),
    legend.box.background = element_rect(color = "black", fill = NA))+
  xlab("Number of cells")+
  ylab("Number of AS genes")

pdf("figure_fghi.merge.20231027.pdf",width = 10,height = 8)
p1/p3|p2/p4
dev.off()


m1<-(p1 / p3) + plot_layout(ncol = 1,heights = c(1,1))
m2<-(p2 / p4) + plot_layout(ncol = 1,heights = c(1,1))

pdf("figure_fghi.merge.v4.pdf",width = 10,height = 8)
wrap_plots(m1,m2,widths = c(1,0.75))
dev.off()












