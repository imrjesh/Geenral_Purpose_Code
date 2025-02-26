setwd("/Users/kumarr9/Desktop/all_sclc_samples_tpm")
library(ggplot2)
#install.packages("readr")
library(readr)
library(data.table)
library(dplyr)
data2 <- list.files(pattern = "*.out")
data3 <- data.frame(data2)
gene_id <- read_tsv("unique_gene.txt",col_names = FALSE)
colnames(gene_id) <- c("Gene_Id")
# n_distinct(gene_id$Gene_Id)
for (i in 1:nrow(data3))
  {
  if(i==1)
    {
    file_rajesh <- read_tsv(paste0("/Users/kumarr9/Desktop/all_sclc_samples_tpm/",data3[i,1]))
    file_rajesh1 <- file_rajesh[,c(1,7)]
    colnames(file_rajesh1) <- c("Gene_Id", paste0(data3[i,1],"@TPM"))
    file_rajesh1 <- file_rajesh1[!duplicated(file_rajesh1$Gene_Id),]
    file_rajesh_2 <- merge(gene_id, file_rajesh1,by=c("Gene_Id"),all = TRUE)
    file_rajesh_2[is.na(file_rajesh_2)]<-0
    write.table(file_rajesh_2, file ="/Users/kumarr9/Downloads/all_rna_sample.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE);
  }
  if(i>1)
    {
    file_rajesh <- read_tsv(paste0("/Users/kumarr9/Desktop/all_sclc_samples_tpm/",data3[i,1]))
    file_rajesh1 <- file_rajesh[,c(1,7)]
    colnames(file_rajesh1) <- c("Gene_Id", paste0(data3[i,1],"@TPM"))
    file_rajesh1 <- file_rajesh1[!duplicated(file_rajesh1$Gene_Id),]
    file_rajesh_2 <- read_tsv("/Users/kumarr9/Downloads/all_rna_sample.tsv");
    file_rajesh_3 <- merge(file_rajesh1, file_rajesh_2,by=c("Gene_Id"),all = TRUE)
    file_rajesh_3[is.na(file_rajesh_3)]<-0
    write.table(file_rajesh_3, file ="/Users/kumarr9/Downloads/all_rna_sample.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE);
  }
}

file_all<-read_tsv("/Users/kumarr9/Downloads/all_rna_sample.tsv");
file_final<-merge(gene_id,file_all,by = c("Gene_Id"))
write.table(file_final, file ="/Users/kumarr9/Downloads/all_rna_sample_final.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE);

                  
                  
              
n_distinct(file_rajesh_3$Gene_Id)
#n_distinct(file_rajesh1$Gene_Id)
#xx <- file_rajesh1[!duplicated(file_rajesh1$Gene_Id),]
