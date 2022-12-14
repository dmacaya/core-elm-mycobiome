## This is the first R script of the analysis. 
## Previously, beginning with the 454 sff files, were run, as standalone programs:
## python process_sff.py, making accesible sff files
## RunTitanium (from AmpliconNoise v.1.29): demultiplex, dechimerize and counting reads reps
## FungalITSextractor
## GramCluster

## GramCluster produces two files, one with the sequences of the clusters (i.e. OTUs): Pooled_GC.fasta
## and another one with the names of the sequences that belong to each cluster: Pooled_GC.text
## The name of the sequences contains info produced by RunTitanium about the times 
## that sequence (or read) was counted in that sample

## This script is to process GramCluster and TunTitanium naming and get total counts.
rm(list=ls())

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/")

db<-read.table("GramCluster_output/Pooled_GC.text",header=F,stringsAsFactors=F)
list_samples<-read.table("Lists/Samples.txt",header=F,stringsAsFactors=F)
list_clusters<-as.numeric(db$V2[!is.na(as.numeric(db$V2))])

count_table<-data.frame(matrix(0,nrow=length(list_clusters),ncol=nrow(list_samples)))
colnames(count_table)<-list_samples$V1
for (i in list_clusters){
  # i<-1
  upper<-match(i,db$V2)+1
  if(i<max(list_clusters)){
    lower<-match(i+1,db$V2)-1
  } else{
    lower<-upper
  }
  temp_db<-db[upper:lower,]
  list_codes<-gsub(">","",temp_db$V2)
  samples<-unlist(strsplit(list_codes,split="_"))[(0:(length(list_codes)-1))*3+1]
  otus<-unlist(strsplit(list_codes,split="_"))[(0:(length(list_codes)-1))*3+2]
  counts<-unlist(strsplit(list_codes,split="_"))[(1:length(list_codes))*3]
  for (j in 1:nrow(list_samples)) {
    # j<-1
    if(sum(samples==list_samples[j,1])>0){
      count_table[i+1,j]<-sum(as.numeric(counts[samples==list_samples[j,1]]))
    } else {
      count_table[i+1,j]<-0
    }
  }
}
count_table<-data.frame(Cluster=paste("Cluster",list_clusters,sep="_"),count_table)

write.table(count_table,"GramCluster_output/Pooled_GC_counts.txt",row.names=F,quote=F,sep="\t")

