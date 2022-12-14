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

sum(counts[,-1])
sum(apply(counts[,-1],1,sum)==1)
sum(apply(counts[,-1],1,sum)==2)
sum(apply(counts[,-1],1,sum)==3)
sum(apply(counts[,-1],1,sum)>5)

table(apply(counts[,-1]>0,1,sum))

sum(!is.na(cluster_info$Species))
sum(!is.na(cluster_info$Genus))
sum(!is.na(cluster_info$Family))
sum(!is.na(cluster_info$Order))
sum(!is.na(cluster_info$Class))
sum(!is.na(cluster_info$Phylum))

min(apply(counts[,-1],2,sum)) 
table(apply(counts[,-1],2,sum))

report<-cluster_info[apply(counts[,-1]>0,1,sum)>7,]
report<-cbind(report,apply(counts[counts$Cluster%in%report$Cluster,-1]>0,1,sum))

report<-cluster_info[apply(counts[,-1]>0,1,sum)>0,]
report<-cbind(report,apply(counts[counts$Cluster%in%report$Cluster,-1]>0,1,sum))

##Collapse by population
wilting<-read.table("Sampling/Wilting.txt",header=T,stringsAsFactors=F)
rivas<-read.table("Sampling/Rivas.txt",header=T,stringsAsFactors=F)
somontes<-read.table("Sampling/Somontes.txt",header=T,stringsAsFactors=F)
somontes<-somontes[1:10,]
burgos<-c("BU22T","BU23T","BU24T")

viv<-apply(counts[,colnames(counts)%in%wilting$Genotipo],1,sum)
riv<-apply(counts[,colnames(counts)%in%rivas$Tree],1,sum)
som<-apply(counts[,colnames(counts)%in%somontes$SAMPLE],1,sum)
bur<-apply(counts[,colnames(counts)%in%burgos],1,sum)

counts_pop<-data.frame(Cluster=cluster_info$Cluster,viv,riv,som,bur)
sum(apply(counts_pop[,-1]>0,1,sum)==1)
sum(apply(counts_pop[,-1]>0,1,sum)==2)
sum(apply(counts_pop[,-1]>0,1,sum)==3)
sum(apply(counts_pop[,-1]>0,1,sum)==4)

sum(viv)
sum(riv)
sum(som)
sum(bur)

sum(viv>0)
sum(riv>0)
sum(som>0)
sum(bur>0)

##Core microbiome - Table 1
VIV_core<-apply(counts[,colnames(counts)%in%wilting$Genotipo]>0,1,sum)>7
SOM_core<-apply(counts[,colnames(counts)%in%somontes$SAMPLE]>0,1,sum)>7
sum(VIV_core)
sum(SOM_core)
sum(VIV_core|SOM_core)

cluster_info_core<-cluster_info[VIV_core|SOM_core,]
counts_core<-counts[VIV_core|SOM_core,]

cluster_info_core_mod<-cluster_info_core[,c(1,3:8)]

output<-data.frame(cluster_info_core_mod,
                   NL=apply(counts[VIV_core|SOM_core,colnames(counts)%in%somontes$SAMPLE]>0,1,sum),
                   NC=apply(counts[VIV_core|SOM_core,colnames(counts)%in%wilting$Genotipo]>0,1,sum),
                   NT=apply(counts[VIV_core|SOM_core,-1]>0,1,sum),
                   Npop=apply(counts_pop[VIV_core|SOM_core,-1]>0,1,sum))
dir.create("Tables")
write.table(output,"Tables/Table_1_Core_microbiome.txt",quote=F,col.names=T,row.names=F)

## Principal Component Analysis
pop_factor<-c(rep("SEBP",10),rep("Rivas",6),rep("Somontes",10),rep("Burgos",3))
population<-data.frame(population=as.factor(pop_factor))
dsdb<-DESeqDataSetFromMatrix(counts_ss[,-1],population,~population)
dsds<-DESeq(dsdb)
vstdb<-varianceStabilizingTransformation(dsdb)

pca<-prcomp(t(assay(vstdb)))
scores<-as.data.frame(pca$x)
pops<-as.factor(c(rep("SEBP",10),rep("Rivas",6),rep("Somontes",10),rep("Burgos",3)))
vars<-pca$sdev^2
perc_vars<-vars/sum(vars)

jpeg("Plots/Fig_2_PCAs.jpg",width=85,height=85,units="mm",res=300)
par(mar=c(5.1,4.1,2.1,2.1),cex=0.6)
plot(scores[,1],scores[,2],
     xlab=paste0("PC1: ",round(perc_vars[1]*100),"% variance"), 
     ylab=paste0("PC2: ",round(perc_vars[2]*100),"% variance"),
     pch=as.numeric(pops))
legend (min(scores[,1]),max(scores[,1])-2,legend=c("Landmark tree","Clonal bank","Rivas population","Burgos population"),
        pch=c(4,3,2,1))
dev.off()

