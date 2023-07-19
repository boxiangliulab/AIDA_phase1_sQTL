###This script was used to generate Figure3B and compare the relative positions of independent cis-sQTLs 
###in the internal or external of introns 

distance<-read.table("C:/users/90410/desktop/aida/replication/cluster_SNP_position_newaida.txt",sep="\t")


frequency<-as.data.frame(matrix(0,110,2))
for(i in 1:nrow(distance)){
  if((distance[i,5]>=distance[i,3])&&(distance[i,5]<=distance[i,4])){
    pro<-(distance[i,5]-distance[i,3])/(distance[i,4]-distance[i,3])
    print(pro)
    print(i)
    if(pro==0){pro<-"0.0"}
    if(distance[i,2]=="+"){
    frequency[51+as.numeric(substr(as.character(pro),3,3)),1]<-frequency[51+as.numeric(substr(as.character(pro),3,3)),1]+1}
    if(distance[i,2]=="-"){
      frequency[60-as.numeric(substr(as.character(pro),3,3)),1]<-frequency[60-as.numeric(substr(as.character(pro),3,3)),1]+1}
  }
  if(distance[i,5]<distance[i,3]){
    pro<-(distance[i,3]-distance[i,5])/500
    pro<-floor(pro)
    if(pro<50){
      if(distance[i,2]=="+"){
    frequency[50-pro,1]<-frequency[50-pro,1]+1}
      if(distance[i,2]=="-"){
        frequency[61+pro,1]<-frequency[61+pro,1]+1}}
  }
  if(distance[i,5]>distance[i,4]){
    pro<-(distance[i,5]-distance[i,4])/500
    pro<-floor(pro)
    if(pro<50){
      if(distance[i,2]=="+"){
    frequency[61+pro,1]<-frequency[61+pro,1]+1}
      if(distance[i,2]=="-"){
        frequency[50-pro,1]<-frequency[50-pro,1]+1}
  }}
}
library(ggplot2)
for(i in 1:50){
  frequency[51-i,2]<-paste(i,"*0.5kb upstream",sep="")
  frequency[60+i,2]<-paste(i,"*0.5kb downstream",sep="")
}
for(i in 1:10){
  frequency[50+i,2]<-paste(i," internal",sep="")
}
colnames(frequency)<-c("Number of sQTL","Position")
frequency$Position<-factor(frequency$Position, levels=frequency$Position)
ggplot(frequency,aes(x=Position,y=`Number of sQTL`))+geom_bar(stat="identity") + theme_classic()+
  theme(title=element_text(size=10),axis.text = element_text(size = 10),axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),legend.title = element_text(size = 7),axis.text.x=element_blank(),axis.ticks.x=element_blank())
