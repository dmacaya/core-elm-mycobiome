rm(list=ls())

library(DESeq2)
library(pheatmap)
library(apeglm)
library(vsn)
library(wesanderson)
library(RColorBrewer)
library(jpeg)
library(shape)
library(hexbin)

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

##Phenotypes
somontes<-read.table("Sampling/Somontes.txt",header=T,stringsAsFactors=F)
somontes<-somontes[1:10,]

## Subsetting to Somontes landmark samples
counts_ss<-counts[,colnames(counts)%in%c("Cluster",somontes$SAMPLE)]
counts_ss<-counts_ss[apply(counts_ss[,-1],1,sum)>0,]

sum(counts_ss[,-1])
sum(apply(counts_ss[,-1],1,sum)==1)
sum(apply(counts_ss[,-1],1,sum)==2)
sum(apply(counts_ss[,-1],1,sum)==3)
sum(apply(counts_ss[,-1],1,sum)>5)

table(apply(counts_ss[,-1]>0,1,sum))

## To create plots
# dir.create("Somontes")
# y<-table(apply(counts_ss[,-1]>0,1,sum))
# save(y,file="Somontes/Somontes_OTU_incidence.Rdata")

cluster_info_ss<-cluster_info[cluster_info$Cluster%in%counts_ss$Cluster,]
sum(!is.na(cluster_info_ss$Family))
sum(!is.na(cluster_info_ss$Order))
sum(!is.na(cluster_info_ss$Class))
sum(!is.na(cluster_info_ss$Phylum))

fam_counts_ss<-fam_counts[,colnames(fam_counts)%in%somontes$SAMPLE]
fam_counts_ss<-fam_counts_ss[apply(fam_counts_ss,1,sum)>0,]
sum(fam_counts_ss)

ord_counts_ss<-ord_counts[,colnames(ord_counts)%in%somontes$SAMPLE]
ord_counts_ss<-ord_counts_ss[apply(ord_counts_ss,1,sum)>0,]
sum(ord_counts_ss)

class_counts_ss<-class_counts[,colnames(class_counts)%in%somontes$SAMPLE]
class_counts_ss<-class_counts_ss[apply(class_counts_ss,1,sum)>0,]
sum(class_counts_ss)

phyl_counts_ss<-phyl_counts[,colnames(phyl_counts)%in%somontes$SAMPLE]
phyl_counts_ss<-phyl_counts_ss[apply(phyl_counts_ss,1,sum)>0,]
sum(phyl_counts_ss)

sum(cluster_info_ss$Phylum=="Ascomycota",na.rm=T)
sum(phyl_counts_ss[1,])/sum(phyl_counts_ss)
sum(cluster_info_ss$Phylum=="Basidiomycota",na.rm=T)
sum(phyl_counts_ss[2,])/sum(phyl_counts_ss)

apply(class_counts_ss,1,sum)
sum(is.na(cluster_info_ss$Class[cluster_info_ss$Phylum=="Ascomycota"]))
sum(is.na(cluster_info_ss$Class[cluster_info_ss$Phylum=="Basidiomycota"]))

temp<-apply(ord_counts_ss,1,sum)
temp[temp>sum(ord_counts_ss)*0.02]

temp<-apply(fam_counts_ss,1,sum)
sort(temp[temp>sum(fam_counts_ss)*0.02],decreasing = T)

apply(counts_ss[,-1],1,sum)

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>7,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>0,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))

