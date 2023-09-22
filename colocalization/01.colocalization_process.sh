###This script was used to perform colocalization analysis to cis-sQTLs of 19 celltypes and 20 GWAS summary statistics.

### First step: select GWAS loci whose p-value < 1e-7 and exclude other GWAS locis whose distance are within 500kb close to the lead GWAS loci.
dir=(
"GBMI_Cell_Genomics_2022_Curated_sumstats_Asthma.txt.gz" \
"HGI_Round7_2022_Curated_sumstats_B2_ALL_eas.txt.gz" \
"Ishigaki_NatGenet_2022_Curated_sumstats_RA_EAS.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_BAS.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_BMI.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_EOS.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_GD.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_Hb.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_Height.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_Ht.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_LYM.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MCH.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MCHC.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MCV.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_MON.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_NEU.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_PLT.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_RBC.txt.gz" \
"Sakaue_Nat_Gen_2021_Curated_sumstats_WBC.txt.gz" \
"Shirai_ARD_2022_Curated_sumstats_AD.txt.gz" \
"Wang_NatCommu_2021_Curated_sumstats_SLE.txt.gz")

for i in ${dir[@]};do
        echo $i
        zcat /ebs1/shared/data/aida/analysis/AIDA_GWAS_sumstats/$i | awk '{ if ($5 < 1e-7) print $0}' \
                > /ebs1/zhangyuntian/sQTL_result/colocalization/test/gwas_overstandard/${i}_gwas_1e_7.txt
done

###R
for(i in 1:length(GWAS)){
        dat<-read.table(paste("../",GWAS[i],"_gwas_1e_7.txt",sep=""),sep="\t")
        pval<-dat[,c(2,5)]
        colnames(pval)<-c("SNP","P")
        pval<-pval[order(pval$P),]
        for(j in 1:(nrow(pval)-1)){
                print(j)
                if(!is.na(pval[j,1])){
                chr<-strsplit(pval[j,1],":")[[1]][1]
                pos<-as.numeric(strsplit(pval[j,1],":")[[1]][2])
                if((chr=="6")&&(pos>=25000000)&&(pos<=35000000))pval[j,1]<-NA
                for(k in (j+1):nrow(pval)){
                        if(!is.na(pval[k,1])){
                        newchr<-strsplit(pval[k,1],":")[[1]][1]
                        newpos<-as.numeric(strsplit(pval[k,1],":")[[1]][2])
                        if(newchr==chr){
                                if((newpos-pos<500000)&&(pos-newpos<500000)){
                                        pval[k,1]<-NA
                                }
                        }}
                }}
        }
        pval<-na.omit(pval)
        write.table(pval,paste(GWAS[i],"_loci",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
        print(GWAS[i])
}


### Second step: identify cis-sQTLs with p-value < 1e-3 && its correponding GWAS loci p-value < 1e-7

celltype=$1

for i in `seq 1 22`;do
        echo $i
        cat /ebs1/zhangyuntian/sQTL_result/nominal_1/${celltype}_${i}_nominals_1.txt | \
                awk '{if($14 < 1e-3) print $0}' > /ebs1/zhangyuntian/sQTL_result/colocalization/test/nominal_overstandard/${celltype}/chr${i}.txt
        awk '{ gsub(/chr/,"", $10); print } ' /ebs1/zhangyuntian/sQTL_result/colocalization/test/nominal_overstandard/${celltype}/chr${i}.txt > \
                /ebs1/zhangyuntian/sQTL_result/colocalization/test/nominal_overstandard/${celltype}/newchr${i}.txt
done

args<-(commandArgs(TRUE))
file_name<-args[1]
celltype<-args[2]
print(celltype)
library(data.table)
GWAS_data<-paste("/ebs1/zhangyuntian/sQTL_result/colocalization/test/gwas_overstandard/",file_name,"_gwas_1e_7.txt",sep="")
print("start read in")
GWAS_pre<-fread(GWAS_data,sep="\t",header=FALSE)
GWAS_pre<-as.data.frame(GWAS_pre)
colnames(GWAS_pre)[1:5]<-c("SNP_hg19","variant_id","beta","se","Pvalue")
GWAS_pre<-GWAS_pre[,c(2,3,4,5)]




##read in nominal_file(from chr1-22)
for(i in 1:22){
  print(paste(i,file_name))
  nominal<-fread(paste("/ebs1/zhangyuntian/sQTL_result/colocalization/test/nominal_overstandard/",celltype,"/newchr",i,".txt",sep=""),sep=" ")
  nominal<-as.data.frame(nominal)
  print("read in nominal")
  if((nrow(nominal)>0)&&(nrow(GWAS_pre)>0)){
    colnames(nominal)[10]<-"variant_id"
    common_snp<-merge(nominal,GWAS_pre,by="variant_id",all=FALSE)
  ##sort common intron and delete 500kb around,exclude chr6:25<A8>C35Mb)
    if(i==6){for(j in 1:nrow(common_snp)){
      seq<-as.numeric(strsplit(common_snp[j,1],":")[[1]][2])
      if((seq>=25000000)&&(seq<=35000000)){common_snp[j,1]<-NA}
    }}
    common_snp<-na.omit(common_snp)
    if(nrow(common_snp)>0){
    common_snp<-common_snp[order(common_snp$Pvalue),]
    for(j in 1:nrow(common_snp)){
      if(is.na(common_snp[j,1])==FALSE){
      seq_j<-as.numeric(strsplit(common_snp[j,1],":")[[1]][2])
      if(nrow(common_snp)>j){
      for(k in (j+1):nrow(common_snp)){
        if(is.na(common_snp[k,1])==FALSE){
        seq_k<-as.numeric(strsplit(common_snp[k,1],":")[[1]][2])
        if((seq_k-seq_j<500000)&&(seq_j-seq_k<500000)){common_snp[k,1]<-NA}
        }
      }}}
    }}
    common_snp<-na.omit(common_snp)
    ##finish delete
    if(nrow(common_snp)>0){
            write.table(common_snp,paste("/ebs1/zhangyuntian/sQTL_result/colocalization/test/inter_loci/",celltype,"/",file_name,i,"loci.txt",sep=""),
                        sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)}
    }
  }

### Third step: build potential datasets with both GWAS loci's information and cis-sQTLs' information
celltype=$1

##build sqtl snp dataset
mkdir inter_loci/${celltype}_nominal_snp
cat /ebs1/zhangyuntian/sQTL_result/colocalization/test/result/${celltype}/new/intron.txt | while read line
do
        echo ${line}
        l=${line%%:*}
        num=${l#*r}
        cat /ebs1/zhangyuntian/sQTL_result/nominal_1/${celltype}_${num}_nominals_1.txt | grep $line | awk '{print $1 "\t" $6 "\t" $10 "\t" $14 "\t" $15}' > \
                /ebs1/zhangyuntian/sQTL_result/colocalization/test/inter_loci/${celltype}_nominal_snp/${line}.txt
done

### Fourth step: perform coloc analysis to the dataset we prepared and obtain the information of PP0-PP4

args<-(commandArgs(TRUE))
file_name<-args[1]
celltype<-args[2]
print(file_name)
library("coloc")
library(dplyr)
library(data.table)

###read in GWAS dataset
GWAS_data<-paste("/ebs1/shared/data/aida/analysis/AIDA_GWAS_sumstats/",file_name,sep="")
print("start read in")
GWAS_pre<-fread(GWAS_data,sep="\t",header=TRUE)
print(paste("read in success",file_name))
GWAS_pre<-as.data.frame(GWAS_pre)
GWAS_pre<-GWAS_pre[,2:5]
colnames(GWAS_pre)<-c("variant_id","beta","se","pval_nominal")
###read in intron need to be dealed
intron_pos<-paste("/ebs1/zhangyuntian/sQTL_result/colocalization/test/inter_loci/intron_each_gwas/",celltype,"_",file_name,".loci",sep="")
loci<-read.table(intron_pos,sep="\t")
###coloc
MAF<-read.table("/ebs1/zhangyuntian/project/aida/genotype/replication/total_MAF.txt",header=TRUE)
colnames(MAF)<-c("variant_id","maf")
print("finished prepare")
if(nrow(loci)>0){
        test_rs<-as.data.frame(matrix(NA,nrow(loci),12))
        colnames(test_rs)<-c("variant_id","intron_id","gene","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","sQTL_pval","gwas_id","gwas_pval")
        for(j in 1:nrow(loci)){
                snp_1mb<-fread(paste("/ebs1/zhangyuntian/sQTL_result/colocalization/test/inter_loci/",celltype,"_nominal_snp/",loci[j,1],".txt",sep=""),sep="\t",header=FALSE)
                colnames(snp_1mb)<-c("gene","intron","variant_id","pval_nominal","beta")
                snp_1mb<-as.data.frame(snp_1mb)
                snp_1mb<-snp_1mb[order(snp_1mb$pval_nominal),]
                test_rs[j,1]<-snp_1mb[1,3]
                test_rs[j,2]<-snp_1mb[1,2]
                test_rs[j,3]<-snp_1mb[1,1]
                test_rs[j,10]<-snp_1mb[1,4]
                for(i in 1:nrow(snp_1mb)){
                  snp_1mb[i,3]<-strsplit(snp_1mb[i,3],"r")[[1]][2]
                }
                input<-merge(snp_1mb,MAF,by="variant_id")
                input<-merge(input,GWAS_pre,by="variant_id",suffixes=c("_sqtl","_gwas"))
                input<-input[order(input$pval_nominal_gwas),]
                se<-input$se_gwas[1]
                test_rs[j,11]<-input$variant_id[1]
                test_rs[j,12]<-input$pval_nominal_gwas[1]
                input$varbeta<-(input$se)^2
                if(nrow(input)>10){
                result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas, type="quant", beta=input$beta_gwas,varbeta=input$varbeta, N=sample_num,snp = input$variant_id),
                          dataset2=list(pvalues=input$pval_nominal_sqtl, type="quant", N=500,snp=input$variant_id), MAF=input$maf)
                test_rs[j,4:9]<-t(as.data.frame(result$summary))[1,1:6]}
        }
        test_rs<-na.omit(test_rs)
        if(nrow(test_rs)>0){
        for(j in 1:nrow(test_rs)){
          if(test_rs[j,9]<=0.75){test_rs[j,1]<-NA}
        }}
        test_rs<-na.omit(test_rs)
  if(nrow(test_rs)>0){write.table(test_rs,paste("coloc_",celltype,file_name,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)}
}
