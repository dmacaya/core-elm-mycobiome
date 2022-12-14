rm(list=ls())

library(DESeq2)
packageVersion("DESeq2")
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

wilting<-read.table("Sampling/Wilting.txt",header=T,stringsAsFactors=F)

########################################SUBSETTING AND BASIC DESCRIPTION#################################
counts_ss<-counts[,colnames(counts)%in%c("Cluster",wilting$Genotipo)]
counts_ss<-counts_ss[apply(counts_ss[,-1],1,sum)>0,]

sum(counts_ss[,-1])
sum(apply(counts_ss[,-1],1,sum)==1)
sum(apply(counts_ss[,-1],1,sum)==2)
sum(apply(counts_ss[,-1],1,sum)==3)
sum(apply(counts_ss[,-1],1,sum)>5)

table(apply(counts_ss[,-1]>0,1,sum))

## To create plots
# dir.create("Vivero")
# y<-table(apply(counts_ss[,-1]>0,1,sum))
# save(y,file="Vivero/Vivero_OTU_incidence.Rdata")

cluster_info_ss<-cluster_info[cluster_info$Cluster%in%counts_ss$Cluster,]
sum(!is.na(cluster_info_ss$Family))
sum(!is.na(cluster_info_ss$Order))
sum(!is.na(cluster_info_ss$Class))
sum(!is.na(cluster_info_ss$Phylum))

fam_counts_ss<-fam_counts[,colnames(fam_counts)%in%wilting$Genotipo]
fam_counts_ss<-fam_counts_ss[apply(fam_counts_ss,1,sum)>0,]
sum(fam_counts_ss)

ord_counts_ss<-ord_counts[,colnames(ord_counts)%in%wilting$Genotipo]
ord_counts_ss<-ord_counts_ss[apply(ord_counts_ss,1,sum)>0,]
sum(ord_counts_ss)

class_counts_ss<-class_counts[,colnames(class_counts)%in%wilting$Genotipo]
class_counts_ss<-class_counts_ss[apply(class_counts_ss,1,sum)>0,]
sum(class_counts_ss)

phyl_counts_ss<-phyl_counts[,colnames(phyl_counts)%in%wilting$Genotipo]
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

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>5,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>3,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))

report<-cluster_info_ss[apply(counts_ss[,-1]>0,1,sum)>0,]
report<-cbind(report,apply(counts_ss[counts_ss$Cluster%in%report$Cluster,-1]>0,1,sum))


##################################ASSOCIATIONS BETWEEN RESISTANCE AND FLORA############################### 

#Phenotypes
divers_trimmed<-divers_trimmed[rownames(divers_trimmed)%in%wilting$Genotipo,]
all(rownames(divers_trimmed)==wilting$Genotipo)
cor.test(divers_trimmed$Shannon,wilting$Wilting)
plot(divers_trimmed$Shannon,wilting$Wilting)
cor.test(divers_trimmed$Simpson,wilting$Wilting)
plot(divers_trimmed$Simpson,wilting$Wilting)
cor.test(divers_trimmed$Rar500,wilting$Wilting)
plot(divers_trimmed$Rar500,wilting$Wilting)

###################OTU_wise######################################
selected_OTUS<-cluster_info_ss[apply(counts_ss[,-1],1,sum)>5,]
counts_selected<-counts_ss[apply(counts_ss[,-1],1,sum)>5,]
wilting<-data.frame(wilting,scl_wilting=scale(wilting$Wilting))
dsdb<-DESeqDataSetFromMatrix(counts_selected[,-1],wilting,~scl_wilting)

##Differential counts
dsds<-DESeq(dsdb,test="Wald")
dsres<-results(dsds)
plotMA(dsres)
dsres_df<-data.frame(dsres)
summary(dsres)
mcols(dsres)$description
dsres_df<-cbind(selected_OTUS,dsres_df)
dsres_df<-dsres_df[order(dsres_df$pvalue),]
write.table(dsres_df,"Stat_analyses/vivero_OTUs_resistance_test.txt",quote=F,row.names = F)

plotMA(dsres)

##Quality of data
plotSparsity(dsds)
ntdb<-normTransform(dsds)
rldb<-rlog(dsds,blind=FALSE)
vstdb<-varianceStabilizingTransformation(dsdb,blind=FALSE)
meanSdPlot(assay(ntdb))
meanSdPlot(assay(rldb))
meanSdPlot(assay(vstdb))
##rldb seems to provide the best variance behavior

dsdf<-estimateSizeFactors(dsdb)
se<-SummarizedExperiment(log2(counts(dsdf, normalized=TRUE) + 1),colData=colData(dsdf))
plotPCA(DESeqTransform(se),intgroup=c("Wilting"))
plotPCA(rldb,intgroup=c("Wilting"))
plotPCA(vstdb,intgroup=c("Wilting"))

##Formal tests
pcaData <- plotPCA(rldb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting) 

pcaData <- plotPCA(vstdb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting)


pheatmap(assay(rldb),cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb),cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE)

hist(rowMeans(counts(dsds,normalized=TRUE)),breaks=100)
select<-rowMeans(counts(dsds,normalized=TRUE))>30&rowMeans(counts(dsds,normalized=TRUE)>0)>0.6
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)

select<-rowMeans(counts(dsds,normalized=TRUE))>20&rowMeans(counts(dsds,normalized=TRUE)>0)>0.6
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)

select<-rowMeans(counts(dsds,normalized=TRUE))>10&rowMeans(counts(dsds,normalized=TRUE)>0)>0.5
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)


###############################Family_wise########################################
fam_counts_ss<-fam_counts[,colnames(fam_counts)%in%wilting$Genotipo]
selected_fam<-fam_counts_ss[apply(fam_counts_ss,1,sum)>5,]
fam_counts_ss<-fam_counts_ss[apply(fam_counts_ss,1,sum)>5,]
dsdb<-DESeqDataSetFromMatrix(fam_counts_ss,wilting,~scl_wilting)

##Differential counts
dsds<-DESeq(dsdb,test="Wald")
dsres<-results(dsds)
plotMA(dsres)
dsres_df<-data.frame(dsres)
summary(dsres)
mcols(dsres)$description
dsres_df<-cbind(selected_fam,dsres_df)
dsres_df<-dsres_df[order(dsres_df$pvalue),]
write.table(dsres_df,"Stat_analyses/vivero_fams_resistance_test.txt",quote=F,row.names = T)

save(list=c("dsres","dsds","dsres_df"),file="Plots/Fig_6_Scatterplots_DED.Rdata")

##Quality of data
plotSparsity(dsds)
ntdb<-normTransform(dsds)
rldb<-rlog(dsds,blind=FALSE)
vstdb<-varianceStabilizingTransformation(dsdb,blind=FALSE)
meanSdPlot(assay(ntdb))
meanSdPlot(assay(rldb))
meanSdPlot(assay(vstdb))
##rldb seems to provide the best variance behavior

dsdf<-estimateSizeFactors(dsdb)
se<-SummarizedExperiment(log2(counts(dsdf, normalized=TRUE) + 1),colData=colData(dsdf))
plotPCA(DESeqTransform(se),intgroup=c("Wilting"))
plotPCA(rldb,intgroup=c("Wilting"))
plotPCA(vstdb,intgroup=c("Wilting"))

pdf("Stat_analyses/PCA_resistance_fam_wise.pdf")
plotPCA(rldb,intgroup=c("Wilting"))
dev.off()

##Formal tests
pcaData <- plotPCA(rldb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting) 

pcaData <- plotPCA(vstdb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting)

pheatmap(assay(rldb),cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb),cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE)

hist(rowMeans(counts(dsds,normalized=TRUE)),breaks=100)
select<-rowMeans(counts(dsds,normalized=TRUE))>30&rowMeans(counts(dsds,normalized=TRUE)>0)>0.6
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)

select<-rowMeans(counts(dsds,normalized=TRUE))>20&rowMeans(counts(dsds,normalized=TRUE)>0)>0.6
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)

select<-rowMeans(counts(dsds,normalized=TRUE))>10&rowMeans(counts(dsds,normalized=TRUE)>0)>0.5
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)


###############################Order_wise########################################
ord_counts_ss<-ord_counts[,colnames(ord_counts)%in%wilting$Genotipo]
selected_ord<-ord_counts_ss[apply(ord_counts_ss,1,sum)>5,]
ord_counts_ss<-ord_counts_ss[apply(ord_counts_ss,1,sum)>5,]
dsdb<-DESeqDataSetFromMatrix(ord_counts_ss,wilting,~scl_wilting)

##Differential counts
dsds<-DESeq(dsdb,test="Wald")
dsres<-results(dsds)
plotMA(dsres)
dsres_df<-data.frame(dsres)
summary(dsres)
mcols(dsres)$description
dsres_df<-cbind(selected_ord,dsres_df)
dsres_df<-dsres_df[order(dsres_df$pvalue),]
write.table(dsres_df,"Stat_analyses/vivero_ords_resistance_test.txt",quote=F,row.names = T)

##Quality of data
plotSparsity(dsds)
ntdb<-normTransform(dsds)
rldb<-rlog(dsds,blind=FALSE)
vstdb<-varianceStabilizingTransformation(dsdb,blind=FALSE)
meanSdPlot(assay(ntdb))
meanSdPlot(assay(rldb))
meanSdPlot(assay(vstdb))
##rldb seems to provide the best variance behavior

dsdf<-estimateSizeFactors(dsdb)
se<-SummarizedExperiment(log2(counts(dsdf, normalized=TRUE) + 1),colData=colData(dsdf))
plotPCA(DESeqTransform(se),intgroup=c("Wilting"))
plotPCA(rldb,intgroup=c("Wilting"))
plotPCA(vstdb,intgroup=c("Wilting"))

##Formal tests
pcaData <- plotPCA(rldb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting) ##Significant

pcaData <- plotPCA(vstdb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting) ##Significant

pheatmap(assay(rldb),cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb),cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE)

hist(rowMeans(counts(dsds,normalized=TRUE)),breaks=100)
select<-rowMeans(counts(dsds,normalized=TRUE))>30&rowMeans(counts(dsds,normalized=TRUE)>0)>0.6
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)

select<-rowMeans(counts(dsds,normalized=TRUE))>10&rowMeans(counts(dsds,normalized=TRUE)>0)>0.5
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
##VAD2 and MDV5 cluster together at order level


###############################Class_wise########################################
class_counts_ss<-class_counts[,colnames(class_counts)%in%wilting$Genotipo]
selected_class<-class_counts_ss[apply(class_counts_ss,1,sum)>5,]
class_counts_ss<-class_counts_ss[apply(class_counts_ss,1,sum)>5,]
dsdb<-DESeqDataSetFromMatrix(class_counts_ss,wilting,~scl_wilting)

##Differential counts
dsds<-DESeq(dsdb,test="Wald")
dsres<-results(dsds)
plotMA(dsres)
dsres_df<-data.frame(dsres)
summary(dsres)
mcols(dsres)$description
dsres_df<-cbind(selected_class,dsres_df)
dsres_df<-dsres_df[order(dsres_df$pvalue),]
write.table(dsres_df,"Stat_analyses/vivero_class_resistance_test.txt",quote=F,row.names = T)

##Quality of data
plotSparsity(dsds)
ntdb<-normTransform(dsds)
rldb<-rlog(dsds,blind=FALSE)
vstdb<-varianceStabilizingTransformation(dsdb,blind=FALSE)
meanSdPlot(assay(ntdb))
meanSdPlot(assay(rldb))
meanSdPlot(assay(vstdb))
##rldb seems to provide the best variance behavior

dsdf<-estimateSizeFactors(dsdb)
se<-SummarizedExperiment(log2(counts(dsdf, normalized=TRUE) + 1),colData=colData(dsdf))
plotPCA(DESeqTransform(se),intgroup=c("Wilting"))
plotPCA(rldb,intgroup=c("Wilting"))
plotPCA(vstdb,intgroup=c("Wilting"))

##Formal tests
pcaData <- plotPCA(rldb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting) ##Significant

pcaData <- plotPCA(vstdb, intgroup=c("Wilting"), returnData=TRUE)
cor.test(pcaData$PC1,pcaData$Wilting)
cor.test(pcaData$PC2,pcaData$Wilting) ##Significant

pheatmap(assay(rldb),cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb),cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE)

hist(rowMeans(counts(dsds,normalized=TRUE)),breaks=100)
select<-rowMeans(counts(dsds,normalized=TRUE))>30&rowMeans(counts(dsds,normalized=TRUE)>0)>0.6
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)

select<-rowMeans(counts(dsds,normalized=TRUE))>10&rowMeans(counts(dsds,normalized=TRUE)>0)>0.5
pheatmap(assay(rldb)[select,],cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE)
pheatmap(assay(rldb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
pheatmap(assay(vstdb)[select,],cluster_rows=F, show_rownames=FALSE,
         cluster_cols=T)
##VAD2 and MDV5 cluster together at class level

