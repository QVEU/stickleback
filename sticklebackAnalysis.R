library(ggplot2)
library(readr)

plotStickle<-function(stickle,name=""){
  pdf(paste("/Volumes/lvd_qve/Sequencing_Data/QVEU_Seq_0009_Minion_InsertionHandleP1_7-20-22/sticklePlots_",name,".pdf",sep = ""))
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

stickleBC1=read_csv("/Volumes/lvd_qve/Sequencing_Data/QVEU_Seq_0009_Minion_InsertionHandleP1_7-20-22/ev71insertion/no_sample/20220718_2031_MC-113212_FAT09008_d7aca3af/fastq_pass/barcode01/merge_stickleback.csv")
plotStickle(stickleBC1,name="BC1")
stickleBC2=read_csv("/Volumes/lvd_qve/Sequencing_Data/QVEU_Seq_0009_Minion_InsertionHandleP1_7-20-22/ev71insertion/no_sample/20220718_2031_MC-113212_FAT09008_d7aca3af/fastq_pass/barcode02/merge_stickleback.csv")
plotStickle(stickleBC2,name="BC2")
ggplot(stickleBC2)+
  stat_summary(aes(insPos_v,...1),fun = "length",cex=0.2)+
  stat_summary(geom = "text",aes(insPos_v,...1,label=insPos_v),fun = "length",hjust=1)+scale_y_log10()

