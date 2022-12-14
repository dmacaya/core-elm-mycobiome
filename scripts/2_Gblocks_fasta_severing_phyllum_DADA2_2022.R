rm(list=ls())

# BiocManager::install("Biostrings")
library(Biostrings)

setwd("D:/Projects/Endophytes/Pyrosequencing/Actions/Ulmus_Paper_0062014_2018/")

cluster_info<-read.table("DADA2/DADA2_Pooled_GC_processed_proper_trimming.txt",header=T,stringsAsFactors=F)

## Loading the output from Gblocks
total_fasta_table<-read.table("Gblocks/Complete_sequences_Clusters_order_curated_2_muscle.fas-gb",
                        header=F,stringsAsFactors=F)
list_clusters<-gsub("[>]","",total_fasta_table[(0:(((nrow(total_fasta_table)-2)/2))*2+1),1])

## Some cluster were removed from the Gblocks run, because of their distant structure
cluster_info<-cluster_info[rownames(cluster_info)%in%list_clusters,]

## Subsetting by taxon and formatting data set
tax_code<-c("taxA","taxB","taxC")
tax_list<-c("Ascomycota","Basidiomycota","Mortierellomycota")
for(i in 1:length(tax_code)){
  taxID<-tax_code[i]
  assign(paste(taxID,"list",sep="_"),
         rownames(cluster_info[cluster_info$Phylum==tax_list[i]&!is.na(cluster_info$Phylum),]))
}
unknown_list<-rownames(cluster_info[is.na(cluster_info$Phylum),])
sum(length(taxA_list)+length(taxB_list)+length(taxC_list)+length(unknown_list))

total_fasta<-readDNAStringSet("Gblocks/Complete_sequences_Clusters_order_curated_2_muscle.fas-gb")

unknown_fasta<-total_fasta[which(list_clusters%in%unknown_list)]
unknown_matrix<-as.matrix(unknown_fasta)
dim(unknown_matrix)

for(i in 1:length(tax_code)){
  taxID<-tax_code[i]
  assign(paste(taxID,"fasta",sep="_"),total_fasta[which(list_clusters%in%get(paste(taxID,"list",sep="_")))])
  assign(paste(taxID,"matrix",sep="_"),as.matrix(get(paste(taxID,"fasta",sep="_"))))
  assign(paste(taxID,"allele_counts",sep="_"),NULL)
  for (i in 1:length(total_fasta[[1]])){
    assign(paste(taxID,"allele_counts",sep="_"),
           rbind(get(paste(taxID,"allele_counts",sep="_")),
                 c(sum(get(paste(taxID,"matrix",sep="_"))[,i]=="A"),
                   sum(get(paste(taxID,"matrix",sep="_"))[,i]=="C"),
                   sum(get(paste(taxID,"matrix",sep="_"))[,i]=="G"),
                   sum(get(paste(taxID,"matrix",sep="_"))[,i]=="T"),
                   sum(get(paste(taxID,"matrix",sep="_"))[,i]=="-"))))
  }
  print(dim(get(paste(taxID,"matrix",sep="_"))))
}

total_matrix<-as.matrix(total_fasta)
apply(total_matrix,2,function(x)sum(x=="-"))


###Extract the segregant SNPS
combined_allele_counts<-NULL
for (i in 1:length(tax_code)){
  taxID<-tax_code[i]
  combined_allele_counts<-cbind(combined_allele_counts,get(paste(taxID,"allele_counts",sep="_")))
}
combined_allele_counts<-combined_allele_counts[
  apply(total_matrix,2,function(x)sum(x=="-"))<nrow(total_matrix)*0.2,]
total_matrix_trimmed<-total_matrix[,apply(total_matrix,2,function(x)sum(x=="-"))<nrow(total_matrix)*0.2]

combined_allele_freq<-NULL
for (i in 0:(length(tax_code)-1)){
  temp<-combined_allele_counts[,(1:4)+i*5]
  combined_allele_freq<-cbind(combined_allele_freq,temp/apply(temp,1,sum))
}

diff<-apply(combined_allele_freq[,(0:(length(tax_code)-1))*4+1],1,function(x)sd(x))+
  apply(combined_allele_freq[,(0:(length(tax_code)-1))*4+2],1,function(x)sd(x))+
  apply(combined_allele_freq[,(0:(length(tax_code)-1))*4+3],1,function(x)sd(x))+
  apply(combined_allele_freq[,(0:(length(tax_code)-1))*4+4],1,function(x)sd(x))

which(diff>0.2)

segregant_alleles<-combined_allele_freq[which(diff>0.2),]

###Generate haplotypes with the segregant alleles
unknown_matrix_trimmed<-unknown_matrix[,apply(total_matrix,2,function(x)sum(x=="-"))<nrow(total_matrix)*0.2]
unknown_matrix_trimmed<-unknown_matrix_trimmed[,which(diff>0.2)]
unknown_haplotypes<-apply(unknown_matrix_trimmed,1,paste,collapse="")
table(unknown_haplotypes)


for(i in 1:length(tax_code)){
  taxID<-tax_code[i]
  assign(paste(taxID,"matrix_trimmed",sep="_"),
         get(paste(taxID,"matrix",sep="_"))[,apply(total_matrix,
                                                   2,function(x)sum(x=="-"))<nrow(total_matrix)*0.2])
  assign(paste(taxID,"matrix_trimmed",sep="_"),
         get(paste(taxID,"matrix_trimmed",sep="_"))[,which(diff>0.2)])
  assign(paste(taxID,"haplotypes",sep="_"),apply(get(paste(taxID,"matrix_trimmed",sep="_")),
                                                 1,paste,collapse=""))
  print(table(get(paste(taxID,"haplotypes",sep="_"))))
}


list_hapotypes<-NULL
for (i in 1:length(tax_code)){
  taxID<-tax_code[i]
  list_hapotypes<-c(list_hapotypes,names(table(get(paste(taxID,"haplotypes",sep="_")))))
}
list_hapotypes<-unique(list_hapotypes)

haplotypes_counts<-NULL
for (j in 1:length(list_hapotypes)){
  temp_vec<-NULL
  for (i in 1:length(tax_code)){
    taxID<-tax_code[i]
    temp_vec<-c(temp_vec,sum(get(paste(taxID,"haplotypes",sep="_"))==list_hapotypes[j]))
  }  
  haplotypes_counts<-rbind(haplotypes_counts,temp_vec)
}
row.names(haplotypes_counts)<-list_hapotypes
haplotypes_counts_trimmed<-haplotypes_counts[apply(haplotypes_counts,1,sum)>1,]


## Creating the assigment rule per haplotype: Threshold set on 0.9, very reliable assigment
assignment_rule<-NULL
for(i in 1:nrow(haplotypes_counts_trimmed)){
  # i<-2
  if(any(haplotypes_counts_trimmed[i,]/sum(haplotypes_counts_trimmed[i,])>0.9)&
     sum(haplotypes_counts_trimmed[i,]>1)){
    for (j in 1:length(tax_code)){
      if(haplotypes_counts_trimmed[i,j]/sum(haplotypes_counts_trimmed[i,])>0.9){
        assignment_rule<-c(assignment_rule,tax_list[j])
      }
    }
  } else {
    assignment_rule<-c(assignment_rule,"Unknown")
  }
}

assignment_rule<-data.frame(haplotypes=rownames(haplotypes_counts_trimmed),
                            assignment_rule,stringsAsFactors = F)

## Assigning taxa to unknown
unknown_assignments<-NULL
for (j in 1:nrow(unknown_matrix_trimmed)){
  # j<-1
  if(unknown_haplotypes[j]%in%assignment_rule$haplotypes){
    unknown_assignments<-c(unknown_assignments,
                           assignment_rule$assignment_rule[assignment_rule$haplotypes==unknown_haplotypes[j]])
  }else{
    unknown_assignments<-c(unknown_assignments,"Unknown")
  }
}
unknown_assignments<-data.frame(ClusterID=rownames(unknown_matrix_trimmed),
                                haplotypes=unknown_haplotypes,assignments=unknown_assignments,
                                stringsAsFactors=F)
sum(unknown_assignments$assignments=="Unknown")

##threshold:0.7; Unknown:1
##threshold:0.8; Unknown:11
##threshold:0.9; Unknown:14


## Filling up taxonomy with the new assignments
cluster_info_post<-read.table("DADA2/DADA2_Pooled_GC_processed_proper_trimming.txt",header=T,stringsAsFactors=F)
sum(cluster_info_post[,"Phylum"]=="Ascomycota",na.rm=T)
sum(cluster_info_post[,"Phylum"]=="Basidiomycota",na.rm=T)
sum(is.na(cluster_info_post[,"Phylum"]))
for(i in 1:nrow(unknown_assignments)){
  # i<-14
  if(unknown_assignments[i,3]=="Unknown"){
    next
  }else {
    cluster_info_post[rownames(cluster_info_post)==unknown_assignments[i,"ClusterID"],"Phylum"]<-unknown_assignments[i,"assignments"]
  }
}
sum(cluster_info_post[,"Phylum"]=="Ascomycota",na.rm=T)
sum(cluster_info_post[,"Phylum"]=="Basidiomycota",na.rm=T)
sum(is.na(cluster_info_post[,"Phylum"]))

write.table(cluster_info_post,"Gblocks/Cluster_taxonomy_phyllum_assign_DADA2.txt",
            quote=F,col.names = T,row.names = T)
write.table(unknown_assignments,"Gblocks/Phyllum_assign_rule_DADA2.txt",
            quote=F,col.names = T,row.names = F)



