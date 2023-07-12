#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Run PCA and Slingshot on B cells 
# author: Yihan Tong
# date: 2023-01-09
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# conda activate r-4.2.3

#######this script helps you to extract dynamic intron junctions across B cell quantiles by ANOVA
### input: all indiduals intron junction ratio file
.libPaths("/ebs1/users/tongyihan/biosoft/miniconda3/envs/sr/bin/Rscript")
library(getopt)
command=matrix(c("file","i",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if( !is.null(args$help) || is.null(args$file))
{ cat(paste(getopt(command,usage=T),"\n"))
  q() }
perind_test=args$file 

print("01 finish read files")

library(tidyverse)
library(stringr)
library(dplyr)
library(reshape2)
#library(onewaytests)

### step 1: calculate all people of 8 quentile

#perind_test<-read.csv("../junction count file/merge_all_sample_test.csv")
#perind_test<-perind_test[,-1]
#for (i in 1:nrow(perind_test)) {
 # b <- as.data.frame(apply(str_split(perind_test[i,2:ncol(perind_test)],"/",simplify = T),2,as.numeric))
 # b <- b[,1]/b[,2]
 # perind_test[i,2:ncol(perind_test)]<-b
#}
#print("02. finish convert fraction to decimal")

#perind_test<-column_to_rownames(perind_test,var = "chrom")

junction<-perind_test$chrom
perind_test_melt<-melt(perind_test,id.vars = "chrom")
perind_test_melt<-dplyr::filter(perind_test_melt,value!="0/0")
perind_test_melt$head<-str_split_fixed(perind_test_melt$value,"\\/",2)[,1]
perind_test_melt$tail<-str_split_fixed(perind_test_melt$value,"\\/",2)[,2]
perind_test_melt$value<-as.numeric(perind_test_melt$head)/as.numeric(perind_test_melt$tail) %>% as.numeric()
perind_test_melt<-perind_test_melt[,c(1:3)]

df<-data.frame()
for (i in junction) {
  df_tmp<-dplyr::filter(perind_test_melt,chrom==i)
  df_tmp$individual<-str_split_fixed(df_tmp$variable,"\\.",3)[,1]
  df_tmp$quantile<-str_split_fixed(df_tmp$variable,"\\.",3)[,2]
  df_tmp$value<-as.numeric(df_tmp$value)
  fit_tmp<-aov(value~quantile,data=df_tmp)
  tabel<-summary(fit_tmp)
  tabel<-as.data.frame(t(unlist(tabel)))
  if (length(grep("Pr",colnames(tabel)))!=0) {
    df_tmp2<-data.frame(junction = i, F_value=tabel$`F value1`,p_value=tabel$`Pr(>F)1`)
  } else {
    df_tmp2<-data.frame(junction = i, F_value=NA,p_value=NA)
  }
  df<-rbind(df,df_tmp2)
}
###remove  introns with "NA" value for both 
delete.na<-function(df,n=0) {
  df[rowSums(is.na(df)) <= n,]
}
df<-delete.na(df)

####extract values from 
index<-which(perind_test$chrom %in% df$junction)
perind_test_extr<-perind_test[index,]
print("03.finish ANOVA")

setwd("02.annova_result")
write.csv(df,"all_sample_ANOVA_value.csv")
write.csv(perind_test_extr,"perind_test_extr.csv")
