###This script was used to plot Figure7A to summary the information of colocalized loci and colocalization events
###upper bar
upper_bar<-read.csv("cell_independent_sqtl_num.csv",header = F)
names(upper_bar)<-c("celltype","number")

upper_bar_sum<-fig7a[,c(1:3)] %>% group_by(cell.type) %>% summarise(sum=sum(number.of.colocalization))


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
                 "CD8\\+_T_GZMB\\+"="GZMB+ CD8+ T",
                 "CD8\\+_T_GZMK\\+"="GZMK+ CD8+ T",
                 "MAIT"="MAIT",
                 "gdT_GZMK-"="GZMK- gdT",
                 "gdT_GZMK\\+"="GZMK+ gdT",
                 "naive_B"="Naive B",
                 "IGHMlo_memory_B"="IGHMlo memory B",
                 "IGHMhi_memory_B"="IGHMhi memory B",
                 "pDC"="pDC",
                 "atypical_B"="Atypical B",
                 "CD4\\+_T" = "CD4+ T") # 22 cell types


for (i in seq_along(label_dict)) {
  upper_bar$cell.type <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],upper_bar$cell.type)
}

upper_bar_sum$cell.type[which(upper_bar_sum$cell.type=="Na茂ve B")] <- "Naive B"
upper_bar_sum$cell.type[which(upper_bar_sum$cell.type=="Na茂ve CD4+ T")] <- "Naive CD4+ T"
upper_bar_sum$cell.type[which(upper_bar_sum$cell.type=="Na茂ve CD8+ T")] <- "Naive CD8+ T"
upper_bar_sum$cell.type[which(upper_bar_sum$cell.type=="GZMB+ CD8+ T ")] <- "GZMB+ CD8+ T"

upper_bar<-merge(upper_bar,upper_bar_sum,by.x="celltype",by.y="cell.type")
upper_bar$percent <- (upper_bar$sum/upper_bar$number)*100
colors2<-c(
  "CD14+ Monocyte"="#fee090", "CD16+ Monocyte"="#fdae61", "CD16+ NK"="#bf812d",
  "CD4+ T"="#4a1486", "cm CD4+ T"="#9e9ac8", "cyt CD4+ T"="#6a51a3", "em CD4+ T"="#807dba",
  "Naive CD4+ T"="#bcbddc", "CD56+ NK"="#8c510a", "GZMB+ CD8+ T"="#3690c0", "GZMK+ CD8+ T"="#74a9cf",
  "Naive CD8+ T"="#a6bddb", "IGHMhi memory B"="#a1d99b", "IGHMlo memory B"="#74c476",
  "MAIT"="#0570b0", "Treg"="#CCCC4D", "Atypical B"="#238b45", "cDC2"="#f46d43",
  "GZMK+ gdT"="#f768a1", "GZMK- gdT"="#dd3497", "Naive B"="#c7e9c0", "pDC"="#d73027")

upper_bar$celltype<-factor(upper_bar$celltype,levels = upper_bar$celltype[rev(order(upper_bar$num_coloc))])
names(upper_bar)<-c("celltype","num_coloc")
ggplot(upper_bar,aes(x=celltype,y=num_coloc,color=celltype,fill=celltype))+
  geom_bar(stat = "identity",width = 0.6)+
  scale_color_manual(values = colors2)+
  scale_fill_manual(values = colors2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +guides(color=FALSE,fill=FALSE)
  

ggsave("figure7A.upper.panel.v2.pdf",width = 8,height = 4)

#middle heatmap
##read in colocalized loci
library(readxl)
colocalized_loci<-read_xlsx("C:/users/90410/desktop/aida/supple_file/table/colocalized_loci_updated_new_7.16.xlsx")
colocalized_loci<-as.data.frame(colocalized_loci)
rownames(colocalized_loci)<-colocalized_loci$...1
colocalized_loci<-colocalized_loci[,-1]
colocalization_event<-read_xlsx("C:/users/90410/desktop/aida/supple_file/table/TableS10_7.17_update.xlsx")
fig7a<-as.data.frame(matrix(NA,20*19,4))
for(i in 1:20){
  for(j in 1:19){
    fig7a[i*19-19+j,1]<-rownames(colocalized_loci)[i]
    fig7a[i*19-19+j,2]<-colnames(colocalized_loci)[j]
    num_1<-which(colocalization_event[,13]==rownames(colocalized_loci)[i])
    num_2<-which(colocalization_event[,14]==colnames(colocalized_loci)[j])
    coloc_num<-length(intersect(num_1,num_2))
    fig7a[i*19-19+j,3]<-coloc_num
    if(colocalized_loci[i,j]=="NA")colocalized_loci[i,j]<-0
    fig7a[i*19-19+j,4]<-as.numeric(colocalized_loci[i,j])/as.numeric(colocalized_loci[i,20])
  }
}
colnames(fig7a)<-c("trait","cell.type","number.of.colocalization","percent.of.colocalized.GWAS.loci")
###read in colocalization events

###
fig7a<-read.csv("AIDA_summary_coloc_results.csv")
fig7a$percent.of.colocalized.GWAS.loci<-100*fig7a$number.of.colocalization/fig7a$number.of.test.loci
fig7a$percent.of.colocalized.GWAS.loci[which(fig7a$percent.of.colocalized.GWAS.loci=="NaN")] <- 0

for(i in 1:20){
  for(j in 1:19){
    fig7a[i*19-19+j,1]<-disease_need[i]
  }
}
fig7a<-fig7a[which(fig7a$trait %in% disease_need),]

for (i in seq_along(label_dict)) {
  fig7a$cell.type <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],fig7a$cell.type)
}
disease_need<-c("Asthma","RA","BAS","BMI","EOS","GD","Hb","Height","Ht","LYM","MCH","MCHC","MCV","MON","NEU","PLT","RBC","WBC","AD","SLE")
fig7a$disease_type<-"a"
fig7a$disease_type[which(fig7a$trait %in% c("AD","Asthma"))] <- "Imflammatory diseases"
fig7a$disease_type[which(fig7a$trait %in% c("GD","RA","SLE"))] <- "Autoimmune diseases"
fig7a$disease_type[which(fig7a$trait %in% c("MON","LYM","EOS","WBC","NEU","BAS"))] <- "Immune related blood traits"
fig7a$disease_type[which(fig7a$trait %in% c("Hb","MCHC","PLT","MCH","MCV","Ht","RBC"))] <- "None immune related diseases"
fig7a$disease_type[which(fig7a$trait %in% c("BMI","Height"))] <- "anthropometric traits"
fig7a$trait<-factor(fig7a$trait,levels = disease_need)

fig7a$cell.type[which(fig7a$cell.type=="Na茂ve B")] <- "Naive B"
fig7a$cell.type[which(fig7a$cell.type=="Na茂ve CD4+ T")] <- "Naive CD4+ T"
fig7a$cell.type[which(fig7a$cell.type=="Na茂ve CD8+ T")] <- "Naive CD8+ T"
fig7a$cell.type[which(fig7a$cell.type=="GZMB+ CD8+ T ")] <- "GZMB+ CD8+ T"

fig7a$cell.type<-factor(fig7a$cell.type,levels = upper_bar$celltype[rev(order(upper_bar$num_coloc))])
fig7a$trait<-factor(fig7a$trait,levels=left_bar$trait[order(left_bar$percent)])
#colors <- colorRampPalette(rev(c("#40004b","#762a83","#9970ab","#c2a5cf","#e7d4e8","#f7f7f7","#d9f0d3","#a6dba0","#5aae61","#1b7837","#00441b")))(500)


colors <- colorRampPalette(rev(c("#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850")))(100)


#strip_col<-c("#d95f02","#7570b3","#e7298a","#e6ab02","#a6761d")
#strip_col<-c("Imflammatory diseases"="#d95f02","Autoimmune diseases"="#7570b3","Immune related blood traits"="#e7298a",
            # "None immune related diseases"="#e6ab02","anthropometric traits"="#a6761d")


ggplot(fig7a, aes(x = cell.type, y = trait)) +
  geom_tile(fill = "white", color = "grey") +
  geom_point(aes(size = percent.of.colocalized.GWAS.loci, fill = number.of.colocalization,color=number.of.colocalization),shape=23,stroke=0.5) +
  labs(x = NULL, y = NULL)+
  scale_fill_gradientn(colours = colors)+
  scale_color_gradientn(colours = colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("figure7A.v4.pdf",width = 8,height = 6)


# left bar
left_bar<-read_xlsx("C:/users/90410/desktop/aida/supple_file/table/percent_colocalized_loci.xlsx")
left_bar<-left_bar[which(left_bar$trait %in% disease_need),]

left_bar$disease_type<-"a"
left_bar$disease_type[which(left_bar$trait %in% c("AD","Asthma"))] <- "Imflammatory diseases"
left_bar$disease_type[which(left_bar$trait %in% c("GD","RA","SLE"))] <- "Autoimmune diseases"
left_bar$disease_type[which(left_bar$trait %in% c("MON","LYM","EOS","WBC","NEU","BAS"))] <- "Immune related blood traits"
left_bar$disease_type[which(left_bar$trait %in% c("Hb","MCHC","PLT","MCH","MCV","Ht","RBC"))] <- "None immune related diseases"
left_bar$disease_type[which(left_bar$trait %in% c("BMI","Height"))] <- "anthropometric traits"
left_bar$trait<-factor(left_bar$trait,levels = left_bar$trait[order(left_bar$percent)])
strip_col<-c("#8da0cb","#e78ac3","#e5c494","#ffd92f","#fc8d62")
strip_col<-c("#fdcdac","#f4cae4","#f1e2cc","#e6f5c9","#fff2ae")


ggplot(left_bar,aes(x=trait,y=percent,color=disease_type,fill=disease_type))+
  geom_bar(stat = "identity",width = 0.6)+
  scale_color_manual(values = strip_col)+
  scale_fill_manual(values = strip_col)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("figure7A.left.panel.pdf",width = 8,height = 3)
