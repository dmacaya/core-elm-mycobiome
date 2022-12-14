rm(list=ls())

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/")

cluster_info<-read.table("Gblocks/Cluster_taxonomy_phyllum_assign_DADA2.txt",header=T,stringsAsFactors=F)
counts<-read.table("GramCluster_output/Pooled_GC_counts_proper_trimming.txt",header=T,stringsAsFactors=F)
counts_trimmed<-counts[counts$Cluster%in%rownames(cluster_info),]
all(cluster_info$Cluster==counts_trimmed$Cluster)
counts<-counts_trimmed
dir.create("Stat_analyses")

##OTU collapse basing on taxonomy
#Family
table(cluster_info$Family)
list_fam<-unique(cluster_info$Family)
list_fam<-list_fam[!is.na(list_fam)]
fam_counts<-NULL
for (i in 1:length(list_fam)){
  # i<-2
  temp_vec<-apply(counts[cluster_info$Family==list_fam[i],-1],2,sum,na.rm=T)
  fam_counts<-rbind(fam_counts,temp_vec)
}
rownames(fam_counts)<-list_fam
write.table(fam_counts,"Stat_analyses/fam_counts_2022.txt",quote=F)
temp_vec<-apply(counts[is.na(cluster_info$Family),-1],2,sum,na.rm=T)
fam_counts<-rbind(fam_counts,temp_vec)
rownames(fam_counts)[nrow(fam_counts)]<-"Unknown"
dim(fam_counts)
write.table(fam_counts,"Stat_analyses/fam_counts_w_unk_2022.txt",quote=F)


#Order
table(cluster_info$Order)
list_ord<-unique(cluster_info$Order)
list_ord<-list_ord[!is.na(list_ord)]
ord_counts<-NULL
for (i in 1:length(list_ord)){
  # i<-2
  temp_vec<-apply(counts[cluster_info$Order==list_ord[i],-1],2,sum,na.rm=T)
  ord_counts<-rbind(ord_counts,temp_vec)
}
rownames(ord_counts)<-list_ord
write.table(ord_counts,"Stat_analyses/ord_counts_2022.txt",quote=F)
temp_vec<-apply(counts[is.na(cluster_info$Order),-1],2,sum,na.rm=T)
ord_counts<-rbind(ord_counts,temp_vec)
rownames(ord_counts)[nrow(ord_counts)]<-"Unknown"
write.table(ord_counts,"Stat_analyses/ord_counts_w_unk_2022.txt",quote=F)

#Class
table(cluster_info$Class)
list_class<-unique(cluster_info$Class)
list_class<-list_class[!is.na(list_class)]
class_counts<-NULL
for (i in 1:length(list_class)){
  # i<-2
  temp_vec<-apply(counts[cluster_info$Class==list_class[i],-1],2,sum,na.rm=T)
  class_counts<-rbind(class_counts,temp_vec)
}
rownames(class_counts)<-list_class
write.table(class_counts,"Stat_analyses/class_counts_2022.txt",quote=F)
temp_vec<-apply(counts[is.na(cluster_info$Class),-1],2,sum,na.rm=T)
class_counts<-rbind(class_counts,temp_vec)
rownames(class_counts)[nrow(class_counts)]<-"Unknown"
write.table(class_counts,"Stat_analyses/class_counts_w_unk_2022.txt",quote=F)

#Phylum
table(cluster_info$Phylum)
list_phyl<-unique(cluster_info$Phylum)
list_phyl<-list_phyl[!is.na(list_phyl)]
phyl_counts<-NULL
for (i in 1:length(list_phyl)){
  # i<-2
  temp_vec<-apply(counts[cluster_info$Phylum==list_phyl[i],-1],2,sum,na.rm=T)
  phyl_counts<-rbind(phyl_counts,temp_vec)
}
rownames(phyl_counts)<-list_phyl
write.table(phyl_counts,"Stat_analyses/phyl_counts_2022.txt",quote=F)
temp_vec<-apply(counts[is.na(cluster_info$Phylum),-1],2,sum,na.rm=T)
phyl_counts<-rbind(phyl_counts,temp_vec)
rownames(phyl_counts)[nrow(phyl_counts)]<-"Unknown"
write.table(phyl_counts,"Stat_analyses/phyl_counts_w_unk_2022.txt",quote=F)
