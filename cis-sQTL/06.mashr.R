###Mashr was used to analyze the proportion of cis-sQTLs that are shared in magnitude
##deal with pvalue(strong set)
##read in intron-SNP pair
gene_SNP<-read.table("C:/users/90410/desktop/aida/mashr/total_pvalue.txt",sep="\t")
num<-which(grepl("chr",gene_SNP[,1],fixed = FALSE)==TRUE)
gene_SNP<-gene_SNP[-num,]
beta_strong<-beta_strong[-num,]
pvalue_strong<-pvalue_strong[-num,]
library(mashr)
celltype_random=c("CD14+_Monocyte","CD16+_Monocyte","CD16+_NK","CD4+_T_cm","CD4+_T_cyt","CD4+_T_em","CD4+_T_naive",
                  "CD56+_NK","CD8+_T_GZMB+","CD8+_T_GZMK+","CD8+_T_naive","IGHMhi_memory_B","IGHMlo_memory_B","MAIT","Treg",
                  "cDC2","gdT_GZMK+","gdT_GZMK-","naive_B")

beta_strong<-read.table("C:/users/90410/desktop/aida/mashr/newbeta_strong.txt",sep="\t")
colnames(beta_strong)<-celltype_random
pvalue_strong<-read.table("C:/users/90410/desktop/aida/mashr/newpvalue_strong.txt",sep="\t")
colnames(pvalue_strong)<-celltype_random


Zscore<-as.data.frame(matrix(0,nrow(pvalue_strong),19))
se_strong<-as.data.frame(matrix(0,nrow(pvalue_strong),19))
for(i in 1:nrow(pvalue_strong)){
  for(j in 1:19){
    if(beta_strong[i,j]>0){Zscore[i,j]<--qnorm(pvalue_strong[i,j]/2)}
    if(beta_strong[i,j]<0){Zscore[i,j]<-qnorm(pvalue_strong[i,j]/2)}
    if(beta_strong[i,j]==0){Zscore[i,j]<-0}
    if(beta_strong[i,j]!=0)se_strong[i,j]<-sqrt(((beta_strong[i,j])^2)/qchisq(pvalue_strong[i,j],1,lower.tail=F))
    if(beta_strong[i,j]==0)se_strong[i,j]<-1000000
  }
  print(i)
}

beta_strong<-as.matrix(beta_strong)
Zscore<-as.matrix(Zscore)
se_strong<-as.matrix(se_strong)
df_list <- list(
              Bhat = beta_strong,
              Shat = se_strong
              )

###read in random dataset
beta<-read.table("C:/users/90410/desktop/aida/mashr/random_beta.txt",sep="\t")
Zscore<-read.table("C:/users/90410/desktop/aida/mashr/random_Zscore.txt",sep="\t")
se<-read.table("C:/users/90410/desktop/aida/mashr/random_se.txt",sep="\t")
df_list_random <- list(
  B = as.matrix(beta),
  Shat = as.matrix(se)
)
data.temp = mash_set_data(df_list_random$Bhat,df_list_random$Shat)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.random = mash_set_data(df_list_random$Bhat,df_list_random$Shat,V=Vhat)
data.strong = mash_set_data(df_list$Bhat,df_list$Shat, V=Vhat)
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

cor<-get_pairwise_sharing(m2) 
write.table(cor,"C:/users/90410/desktop/aida/supple_file/corr_aida.txt",row.names=TRUE,col.names=TRUE,quote=FALSE)


###Number of tissues shared by sign and magnitude

beta_strong<-read.table("C:/users/90410/desktop/aida/mashr/newbeta_strong.txt",sep="\t")
pvalue_strong<-read.table("C:/users/90410/desktop/aida/mashr/newpvalue_strong.txt",sep="\t")
beta_strong<-as.matrix(beta_strong)
pvalue_strong<-as.matrix(pvalue_strong)
zscore<--qnorm(pvalue_strong/2)*beta_strong/abs(beta_strong)
colnames(beta_strong)<-celltype_random
sign.norm = function(effectsize) {
  t(apply(effectsize,1,function(x){
    x/sign(x[which.max(abs(x))])
  }))}

sign.tissue.func = function(normdat){
  apply(normdat,1,function(x){
    sum(x<0)})}
sign.func <- function (normeffectsize)
  apply(normeffectsize,1,function(x)(sum(x>0)))
het.norm = function(effectsize) {
  t(apply(effectsize,1,function(x){
    x/x[which.max(abs(x))]
  }))
}

het.func = function (normdat, threshold) {
  apply((normdat),1,function(x){sum(x > threshold)})
}

lfsr<-read.table("C:/users/90410/desktop/aida/mashr/LFSR_4.23.txt",sep="\t",header=TRUE)

posterior.mean<-mash_compute_posterior_matrices(m2,data)$PosteriorMean

pm.mash.beta<-posterior.mean*beta_strong/zscore

###plot the two together
magnitude<-het.func(het.norm(effectsize=beta_strong[nsig>0,]),threshold=0.5)
sign<-sign.func(het.norm(effectsize=beta_strong[nsig>0,]))

toge<-as.data.frame(matrix(NA,19,3))
toge$V1<-c(1:19)
for(i in 1:19){
  num_1<-which(magnitude[]==i)
  num_2<-which(sign[]==i)
  toge[i,2]<-length(num_1)/length(magnitude)
  toge[i,3]<-length(num_2)/length(sign)
}
colnames(toge)<-c("Number of cell types","magnitude","sign")

##画分组柱状图
bar_plot<-data.frame(`Number of cell types`=rep(c(1:19),2),Fraction=append(toge$magnitude,toge$sign),type=c(rep("Magnitude",19),rep("Sign",19)))
bar_plot$Number.of.cell.types<-factor(bar_plot$Number.of.cell.types,levels=c(1:19))
ggplot(bar_plot,aes(x=Number.of.cell.types,y=Fraction,fill=type))+
  geom_bar(stat="identity",position="dodge",width=0.8,color="black")+
  xlab("Number of cell types")+ylab("Fraction")+theme_classic()+theme(legend.title=element_blank(),legend.position = "top")
