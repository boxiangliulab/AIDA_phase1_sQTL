#   ____________________________________________________________________________
#   Script information                                                      ####

# title: Calculate intron type proportion in Snaptrons
# author: Yihan Tong
# date: 2023-03-10
# description: None

#   ____________________________________________________________________________
#   HPC details                                                             ####

# conda activate r-4.2.3

#' this script is to calculate intron type proportion in Snaptrons
#' input: junctions file downloaded from Snaptrons srav2 and all test intron junction

library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)

snaptrons<-fread("junctions",sep="\t",header=F)
print("01.finish read snaptrons")

ref_junction<-data.frame(chr=snaptrons$V2,start=snaptrons$V3,end=snaptrons$V4,strand=snaptrons$V6,sample=snaptrons$V12,sample_count=snaptrons$V13)
print("02.finish wite ref_junction_useful_info")

junctions<-read.csv("aida_perind.counts.gz.qqnorm_total_no_chrX_chr_Y_only_pc_lncRNA_only_position_with_strand_1_based_v2.txt",sep=" ")
junctions$strand<-str_split_fixed(junctions$ID,"_",3)[,3]

for (num in c(1,10,100,300,500,1000)) {
  ref_junction<-data.frame(chr=snaptrons$V2,start=snaptrons$V3,end=snaptrons$V4,strand=snaptrons$V6,sample=snaptrons$V12,sample_count=snaptrons$V13)
  ref_junction<-dplyr::filter(ref_junction,sample_count>=num)
  setwd("~result/")
  dir.create(as.character(num))
  setwd(paste0("~result/",num))
  ####for + strand
  ref_junction <- dplyr::filter(ref_junction,chr != "chrX" & chr != "chrY")
  data<-dplyr::filter(junctions,strand=="+")
  ref_junction<-dplyr::filter(ref_junction,strand=="+")
  
  merge_result<-data.frame()
  names(data)<-c("chr","start","end","clu","strand")
  junction_num<-as.numeric(nrow(data))
  both_present_plus<-data.frame()
  both_absent_plus<-data.frame()
  donor_absent_plus<-data.frame()
  acceptor_absent_plus<-data.frame()
  for (chromosome in unique(data$chr)) {
    data_chr<-dplyr::filter(data,chr==chromosome)
    ref_chr<-dplyr::filter(ref_junction,chr==chromosome)
    for (i in 1:nrow(data_chr)) {
      if ((data_chr$start)[i] %in% ref_chr$start) {
        ref_index<-which(ref_chr$start==(data_chr$start)[i])
        if ((data_chr$end)[i] %in% ref_chr$end[ref_index]) {
          both_present_tmp<-data_chr[i,]#both present
          both_present_plus<-rbind(both_present_plus,both_present_tmp)
        } else {
          acceptor_absent_tmp<-data_chr[i,]#acceptor_absent
          acceptor_absent_plus<-rbind(acceptor_absent_tmp,acceptor_absent_plus)
        }
      } else {
        if ((data_chr$end)[i] %in% ref_chr$end) {
          donor_present_tmp<-data_chr[i,]#donor absent
          donor_absent_plus<-rbind(donor_absent_plus,donor_present_tmp)
        } else {
          both_absent_tmp<-data_chr[i,]
          both_absent_plus<-rbind(both_absent_plus,both_absent_tmp)#both absent
        }
      }
    }
  }
  
  # for - strand
  ref_junction<-data.frame(chr=snaptrons$V2,start=snaptrons$V3,end=snaptrons$V4,strand=snaptrons$V6,sample=snaptrons$V12,sample_count=snaptrons$V13)
  ref_junction <- dplyr::filter(ref_junction,chr != "chrX" & chr != "chrY")
  #ref_junction<-dplyr::filter(ref_junction,V4=="-")
  data<-dplyr::filter(junctions,strand=="-")
  ref_junction<-dplyr::filter(ref_junction,strand=="-")
  merge_result<-data.frame()
  names(data)<-c("chr","start","end","clu","strand")
  junction_num<-as.numeric(nrow(data))
  both_present_minus<-data.frame()
  both_absent_minus<-data.frame()
  donor_absent_minus<-data.frame()
  acceptor_absent_minus<-data.frame()
  for (chromosome in unique(data$chr)) {
    data_chr<-dplyr::filter(data,chr==chromosome)
    ref_chr<-dplyr::filter(ref_junction,chr==chromosome)
    for (i in 1:nrow(data_chr)) {
      if ((data_chr$end)[i] %in% ref_chr$end) {
        ref_index<-which(ref_chr$end==(data_chr$end)[i])
        if ((data_chr$start)[i] %in% ref_chr$start[ref_index]) {
          both_present_tmp<-data_chr[i,]#both present
          both_present_minus<-rbind(both_present_minus,both_present_tmp)
        } else {
          acceptor_absent_tmp<-data_chr[i,]#acceptor_absent
          acceptor_absent_minus<-rbind(acceptor_absent_tmp,acceptor_absent_minus)
        }
      } else {
        if ((data_chr$start)[i] %in% ref_chr$start) {
          donor_present_tmp<-data_chr[i,]#donor absent
          donor_absent_minus<-rbind(donor_absent_minus,donor_present_tmp)
        } else {
          both_absent_tmp<-data_chr[i,]
          both_absent_minus<-rbind(both_absent_minus,both_absent_tmp)#both absent
        }
      }
    }
  }
  
  acceptor_absent<-rbind(acceptor_absent_minus,acceptor_absent_plus)
  donor_absent<-rbind(donor_absent_minus,donor_absent_plus)
  both_absent<-rbind(both_absent_minus,both_absent_plus)
  both_present<-rbind(both_present_minus,both_present_plus)
  
  
  
  write.csv(acceptor_absent,"./acceptor_absent.csv",row.names = F)
  write.csv(donor_absent,"./donor_absent.csv",row.names = F)
  write.csv(both_absent,"./both_absent.csv",row.names = F)
  write.csv(both_present,"./both_present.csv",row.names = F)
  merge_result_tmp<-data.frame(donor_absent=(nrow(donor_absent_plus)+nrow(donor_absent_minus))/nrow(junctions),
                               acceptor_absent=(nrow(acceptor_absent_minus)+nrow(acceptor_absent_plus))/nrow(junctions),
                               both_absent=(nrow(both_absent_minus)+nrow(both_absent_plus))/nrow(junctions),
                               both_present=(nrow(both_present_minus)+nrow(both_present_plus))/nrow(junctions))
  merge_result_tmp$sum_result=sum(merge_result_tmp[1,c(1:4)])
  merge_result_tmp$junction_number<-nrow(junctions)
  merge_result<-rbind(merge_result,merge_result_tmp)
  write.csv(merge_result,"./merge_result.csv",row.names = F)
  
}