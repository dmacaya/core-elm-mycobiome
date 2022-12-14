rm(list=ls())

library(dada2); packageVersion("dada2")

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/") 

sequences<-read.table("GramCluster_output/Pooled_GC.fasta",header=F,stringsAsFactors=F,sep="\t")
temp_seq_1<-sequences[1+(0:((nrow(sequences)-1)/2))*2,]
temp2_seq_1<-unlist(strsplit(temp_seq_1,split=", "))
seq_1<-gsub("[>]","",temp2_seq_1[1+(0:((nrow(sequences)-1)/2))*2])
seq_1<-gsub(" ","_",seq_1)
seq_2<-sequences[2+(0:((nrow(sequences)-1)/2))*2,]

counts<-read.table("GramCluster_output/Pooled_GC_counts.txt",header=T,stringsAsFactors=F)
all(counts$Cluster==seq_1)
sum(counts[,-c(1,28,29)])

seqtab<-t(counts[,-1])
colnames(seqtab)<-seq_2

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

which(!colnames(seqtab)%in%colnames(seqtab.nochim))

## Download UNITE database from their website and place it in folder blastn_output
unite.ref <- "blastn_output/sh_general_release_dynamic_all_16.10.2022.fasta"  
taxa <- assignTaxonomy(seqtab, unite.ref, multithread = TRUE, tryRC = TRUE)
taxa_df<-data.frame(taxa)
all(rownames(taxa_df)==colnames(seqtab))

taxa_print<-taxa_df
rownames(taxa_print)<-counts$Cluster

taxa_print$Kingdom<-gsub("k__","",taxa_print$Kingdom)
taxa_print$Phylum<-gsub("p__","",taxa_print$Phylum)
taxa_print$Class<-gsub("c__","",taxa_print$Class)
taxa_print$Order<-gsub("o__","",taxa_print$Order)
taxa_print$Family<-gsub("f__","",taxa_print$Family)
taxa_print$Genus<-gsub("g__","",taxa_print$Genus)
taxa_print$Species<-gsub("s__","",taxa_print$Species)

cluster_to_keep<-rownames(taxa_print)[taxa_print$Kingdom=="Fungi"&!is.na(taxa_print$Phylum)]
cluster_to_remove<-rownames(taxa_print)[taxa_print$Kingdom!="Fungi"|is.na(taxa_print$Phylum)]

taxa_print_trimmed<-taxa_print[rownames(taxa_print)%in%cluster_to_keep,]

dir.create("DADA2")
write.table(taxa,"DADA2/DADA2_Pooled_GC_raw_output.txt",
            quote=F,col.names=T,row.names=T)
write.table(taxa_print,"DADA2/DADA2_Pooled_GC_processed_output.txt",
            quote=F,col.names=T,row.names=T)
write.table(taxa_print_trimmed,"DADA2/DADA2_Pooled_GC_trimmed_output.txt",
            quote=F,col.names=T,row.names=T)
write.table(seqtab,"DADA2/DADA2_Pooled_GC_seqtab.txt",
            quote=F,col.names=T,row.names=T)
write.table(seqtab.nochim,"DADA2/DADA2_Pooled_GC_seqtab_nochim.txt",
            quote=F,col.names=T,row.names=T)
