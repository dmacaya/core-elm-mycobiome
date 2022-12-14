rm(list=ls())

library(DESeq2)
library(pheatmap)
library(apeglm)
library(vsn)
library(wesanderson)
library(RColorBrewer)
library(jpeg)
library(shape)

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/")

cluster_info<-read.table("Gblocks/Cluster_taxonomy_phyllum_assign_DADA2.txt",header=T,stringsAsFactors=F)
cluster_info<-data.frame(Cluster=rownames(cluster_info),cluster_info)
counts<-read.table("GramCluster_output/Pooled_GC_counts_proper_trimming.txt",header=T,stringsAsFactors=F)
counts_trimmed<-counts[counts$Cluster%in%rownames(cluster_info),]
all(cluster_info$Cluster==counts_trimmed$Cluster)
counts<-counts_trimmed

fam_counts<-read.table("Stat_analyses/fam_counts_2022.txt",header=T ,stringsAsFactors=F)
ord_counts<-read.table("Stat_analyses/ord_counts_2022.txt",header=T ,stringsAsFactors=F)
class_counts<-read.table("Stat_analyses/class_counts_2022.txt",header=T ,stringsAsFactors=F)
phyl_counts<-read.table("Stat_analyses/phyl_counts_2022.txt",header=T ,stringsAsFactors=F)
divers_trimmed<-read.table("Diversity/Diversity_trimmed_abs.txt",header=T ,stringsAsFactors=F)

divers_trimmed<-divers_trimmed[rownames(divers_trimmed)%in%rivas$Tree,]

## Phenotypes and subset
rivas<-read.table("Sampling/Rivas.txt",header=T,stringsAsFactors=F)

counts_ss<-counts[,colnames(counts)%in%c("Cluster",rivas$Tree)]
counts_ss<-counts_ss[apply(counts_ss[,-1],1,sum)>0,]

sum(counts_ss[,-1])
sum(apply(counts_ss[,-1],1,sum)==1)
sum(apply(counts_ss[,-1],1,sum)==2)
sum(apply(counts_ss[,-1],1,sum)==3)
sum(apply(counts_ss[,-1],1,sum)>5)

table(apply(counts_ss[,-1]>0,1,sum))

temp<-counts_ss[apply(counts_ss[,-1]>0,1,sum)==5,]

## To create plots
# dir.create("Rivas")
# y<-table(apply(counts_ss[,-1]>0,1,sum))
# save(y,file="Rivas/Rivas_OTU_incidence.Rdata")

cluster_info_ss<-cluster_info[cluster_info$Cluster%in%counts_ss$Cluster,]
sum(!is.na(cluster_info_ss$Species))
sum(!is.na(cluster_info_ss$Genus))
sum(!is.na(cluster_info_ss$Family))
sum(!is.na(cluster_info_ss$Order))
sum(!is.na(cluster_info_ss$Class))
sum(!is.na(cluster_info_ss$Phylum))

fam_counts_ss<-fam_counts[,colnames(fam_counts)%in%rivas$Tree]
fam_counts_ss<-fam_counts_ss[apply(fam_counts_ss,1,sum)>0,]
sum(fam_counts_ss)

ord_counts_ss<-ord_counts[,colnames(ord_counts)%in%rivas$Tree]
ord_counts_ss<-ord_counts_ss[apply(ord_counts_ss,1,sum)>0,]
sum(ord_counts_ss)

class_counts_ss<-class_counts[,colnames(class_counts)%in%rivas$Tree]
class_counts_ss<-class_counts_ss[apply(class_counts_ss,1,sum)>0,]
sum(class_counts_ss)

phyl_counts_ss<-phyl_counts[,colnames(phyl_counts)%in%rivas$Tree]
phyl_counts_ss<-phyl_counts_ss[apply(phyl_counts_ss,1,sum)>0,]
sum(phyl_counts_ss)

min(apply(counts_ss[,-1],2,sum)) 
table(apply(counts_ss[,-1],2,sum))

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>4,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>0,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))

