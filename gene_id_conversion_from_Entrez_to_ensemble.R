#To Convert Entrez gene ids to Ensembl gene ids

#Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

#Selecting a BioMart database and dataset
library("biomaRt")
listMarts()
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#How to build a biomaRt query
#The getBM() function has three arguments that need to be introduced: filters, attributes and values.
filters = listFilters(ensembl)
filters[1:5,] #to show a few list of filters

#values : input data file
entrezgene = c("3098","728642") 

#attributes : type of output 

#getBM() function is the main query function in biomaRt.
genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id"), values=entrezgene, mart=ensembl)                                                                                                                 
print(genes)

write.table(genes, paste(getwd(), "/output.txt", sep=""), quote=FALSE, append=FALSE)
