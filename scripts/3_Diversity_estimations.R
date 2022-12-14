rm(list=ls())

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/")
require(vegan)
packageVersion("vegan")

cluster_info_post<-read.table("Gblocks/Cluster_taxonomy_phyllum_assign_DADA2.txt",header=T,stringsAsFactors=F)
counts<-read.table("GramCluster_output/Pooled_GC_counts_proper_trimming.txt",header=T,stringsAsFactors=F)
counts_trimmed<-counts[counts$Cluster%in%rownames(cluster_info_post),]

dat.trimmed <- data.frame(t(counts_trimmed[,-1]))
colnames(dat.trimmed)<-counts_trimmed$Cluster

shannon.trimmed <- diversity (dat.trimmed, index = "shannon")
shannon.trimmed <- as.data.frame(shannon.trimmed)
names(shannon.trimmed) <- c("Shannon")
simpson.trimmed <- diversity (dat.trimmed, index = "simpson")
simpson.trimmed <- as.data.frame(simpson.trimmed)
names(simpson.trimmed) <- c("Simpson")
rar.trimmed <- rarefy (dat.trimmed, 500,se=T)
rar.trimmed <- as.data.frame(t(rar.trimmed))
names(rar.trimmed) <- c("Rar500","SE")
total.trimmed <- cbind (shannon.trimmed, simpson.trimmed, rar.trimmed)
dir.create("Diversity")
write.table(total.trimmed, "Diversity/Diversity_trimmed_abs.txt", sep="\t")

x<-rarecurve(dat.trimmed, 50)
limit<-max(unlist(lapply(x,length)))-1
x_axis<-attributes(x[[which(lapply(x,length)==max(unlist(lapply(x,length))))]])$
  Subsample[-max(unlist(lapply(x,length)))]
output<-NULL
for (i in 1:nrow(dat.trimmed)){
  # i<-1
  temp_vec<-x[[i]][1:(length(x[[i]])-1)]
  temp_vec_fill<-c(temp_vec,rep(NA,limit-length(temp_vec)))
  output<-cbind(output,temp_vec_fill)
}
colnames(output)<-row.names(dat.trimmed)

dir.create("Plots")
save(list=c("x","x_axis","output"),file="Plots/Fig_2_Rarefaction_trimmed.Rdata")
