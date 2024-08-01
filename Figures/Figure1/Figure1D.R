#' this script allows you to visualize the result of intron type of snaptrons and gencode
gen_merge<-read.csv("merge_result.csv")
gen_merge$type<-"Gencode"
snap_merge_1<-read.csv("merge_result.csv")
snap_merge_1$type<-"Snaptrons(N=1)"
snap_merge_10<-read.csv("merge_result.csv")
snap_merge_10$type<-"Snaptrons(N=10)"
snap_merge_100<-read.csv("merge_result.csv")
snap_merge_100$type<-"Snaptrons(N=100)"
snap_merge_300<-read.csv("merge_result.csv")
snap_merge_300$type<-"Snaptrons(N=300)"
snap_merge_500<-read.csv("merge_result.csv")
snap_merge_500$type<-"Snaptrons(N=500)"

snap_merge_1000<-read.csv("merge_result.csv")
snap_merge_1000$type<-"Snaptrons(N=1000)"

merge_result<-rbind(gen_merge,snap_merge_1,snap_merge_100,snap_merge_10,snap_merge_300,snap_merge_500,snap_merge_1000)
merge_result<-merge_result[,-c(5,6)]
merge_result<-merge_result[order(-merge_result$both_present),]
merge_result$one_end_mising<-merge_result$donor_absent+merge_result$acceptor_absent # merge acceptor missing and donor missing
merge_result<-merge_result[,-c(1,2)]
#merge_result$both_present<-factor(merge_result$both_present,levels =merge_result$both_present )

##add long read data
pacbi<-read.csv("pacbio4 individuals.csv")
pacbio<-data.frame(both_absent=pacbi$both_absent,both_present=pacbi$both_present,type=paste0("pacbio_N=",pacbi$individual_num),one_end_mising=pacbi$donor_absent+pacbi$acceptor_absent)
#merge long read data with gencode and snaptron

merge_result<-rbind(merge_result,pacbio)

library(reshape2)
result<-melt(merge_result)

#box1 is the result with gencode and pacbio
box1<-result[c(grep("Gencode",result$type),grep("pacbio",result$type)),]
box1 <- dplyr::filter(box1,variable=="both_present")
none<-data.frame(type="none",variable="both_present",value=1)
box1<-rbind(box1,none)
box1$difference<-c(box1$value[1],diff(box1$value))
box1$type<-factor(box1$type,levels = rev(c("Gencode","pacbio_N=1","pacbio_N=2","pacbio_N=3","pacbio_N=4","none")))


pdf("Figure1D_proportion_1.pdf",height = 1,width = 4)
ggplot(box1,aes(x=difference,y=variable,fill=type))+
  geom_bar(stat = "identity")+
  labs(title = "Total intron junctions: 50102")+
  scale_fill_manual(values = c(none="#f781bf","pacbio_N=2"="#a65628","pacbio_N=1"="#377eb8","pacbio_N=4"="#984ea3","pacbio_N=3"="#ffff33","Gencode"="#e41a1c"))+
  ylab("GENCODE & PacBio")+
  xlab("Proportion")+
  theme_classic()
dev.off()

#box2 is the result with snaptron
box2<-result[c(grep("Snaptrons",result$type)),]
box2 <- dplyr::filter(box2,variable=="both_present")
box2<-box2[c(6,5,4,3,2,1),]
none<-data.frame(type="none",variable="both_present",value=1)
box2<-rbind(box2,none)
box2$relative<-box2$value-0.9
box2$value<-box2$relative/0.1
box2$difference<-c(box2$value[1],diff(box2$value))

box2$type<-factor(box2$type,levels = rev(c("Snaptrons(N=1000)","Snaptrons(N=500)","Snaptrons(N=300)","Snaptrons(N=100)","Snaptrons(N=10)","Snaptrons(N=1)","none")))

pdf("Figure1D_proportion_2.pdf",height = 1,width = 4)
ggplot(box2,aes(x=difference,y=variable,fill=type))+
  geom_bar(stat = "identity")+
  labs(title = "Total intron junctions: 50102")+
  scale_fill_manual(values = c(none="#d9d9d9","Snaptrons(N=1)"="#8dd3c7","Snaptrons(N=10)"="#bebada","Snaptrons(N=100)"="#ffffb3","Snaptrons(N=300)"="#fccde5","Snaptrons(N=500)"="#80b1d3","Snaptrons(N=1000)"="#fb8072"))+
  ylab("Snaptron")+
  xlab("Proportion")+
  theme_classic()
dev.off()
