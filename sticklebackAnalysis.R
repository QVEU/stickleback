library(ggplot2)
library(readr)
inputdir="/Volumes/LVD_QVE/Sequencing_Data/QVEU_Seq_0018_Short_Long_Insertions/shortinsertions/no_sample/20220913_2012_MC-113212_FAT29354_04a6d4cf/fastq_pass/barcode24"
plotStickle<-function(inputdir,name=""){
  stickle=read_csv(paste(inputdir,"/merge_stickleback.csv", sep = ""))
  pdf(paste(inputdir,"/sticklePlots_",name,".pdf",sep = ""))
  print(stickle$minD_v)
  plot(ggplot(stickle)+geom_histogram(aes(minD)))
  plot(ggplot(stickle)+geom_histogram(aes(insPos_v,fill=factor(minD)),position="stack",bins =100))
  plot(ggplot(stickle)+geom_histogram(aes(insPos_v,fill=factor(minD)),position="stack",binwidth=1))
  plot(ggplot(stickle)+
         stat_summary(aes(insPos_v,...1),fun = "length",cex=0.2)+
         stat_summary(geom = "text",aes(insPos_v,...1,label=insPos_v),fun = "length",hjust=1))
  try(heatmap(scale = "none",Rowv = NA,Colv = NA,table(stickle$minD,stickle$minD_v)))
  tableMP<-table(stickle$minPos)
  tableMP<-data.frame(tableMP)
  plot(ggplot(tableMP)+geom_line(aes(rank(-Freq),Freq/sum(Freq)))+scale_y_log10()+ylab("Insertion Frequency")+xlab("Position Rank"))
  dev.off()
}

plotStickle(inputdir,name="Insertions")

