rm(list=ls())

library(dada2); packageVersion("dada2")

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/")

taxa_print<-read.table("DADA2/DADA2_Pooled_GC_processed_output.txt",header=T,stringsAsFactors=F)
taxa_print_trimmed<-read.table("DADA2/DADA2_Pooled_GC_trimmed_output.txt",header=T,stringsAsFactors=F)
seqtab<-read.table("DADA2/DADA2_Pooled_GC_seqtab.txt",header=T,stringsAsFactors=F)
seqtab_nochim<-read.table("DADA2/DADA2_Pooled_GC_seqtab_nochim.txt",header=T,stringsAsFactors=F)

cluster_info<-taxa_print

counts<-read.table("GramCluster_output/Pooled_GC_counts.txt",header=T,stringsAsFactors=F)

counts_ss<-counts[,!colnames(counts)%in%c("H1N","H1S")] ## These samples were just internal checks and never used in the analyses
index_som_extra<-which(apply(counts_ss[,-1],1,sum)==0)

index_chimera<-which(!colnames(seqtab)%in%colnames(seqtab_nochim))

table(cluster_info$Phylum)
table(cluster_info$Kingdom)
index_rare<-which(cluster_info$Kingdom%in%c("Viridiplantae","Rhizaria","Metazoa"))

table(cluster_info$Class)

index_singletons<-which(apply(counts_ss[,-1],1,sum)==1)
index_doubletons<-which(apply(counts_ss[,-1],1,sum)==2)
index_tripletons<-which(apply(counts_ss[,-1],1,sum)==3)

sum(!index_singletons%in%index_som_extra)
sum(apply(counts_ss[,-1],1,sum)>5)

index_DADA_trimmed<-which(!rownames(taxa_print)%in%rownames(taxa_print_trimmed))

joint_index<-unique(c(index_chimera,index_rare,index_singletons,index_doubletons,index_som_extra))

cluster_info_trimmed<-cluster_info[-joint_index,]

counts_trimmed<-counts_ss[counts_ss$Cluster%in%rownames(cluster_info_trimmed),]

sum(counts[,-1])
sum(counts_trimmed[,-1])

sum(counts_trimmed[,-1])/sum(counts[,-1])

write.table(cluster_info_trimmed,"DADA2/DADA2_Pooled_GC_processed_proper_trimming.txt",
            quote=F,col.names = T,row.names = T)
write.table(counts_trimmed,"GramCluster_output/Pooled_GC_counts_proper_trimming.txt",
            quote=F,col.names = T,row.names = F)

