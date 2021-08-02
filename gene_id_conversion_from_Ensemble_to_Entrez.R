library(biomaRt)

#example

#input
ensembl.genes = c("ENSG00000175899") #a vector of Ensembl gene IDs

#change
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id"),
    values=ensembl.genes,
    mart=mart)



#real data

#input
input_data = read.table("C:/R/default_working_directory/Secretome/Human Protein Atlas data/sa_location_Secreted.tsv", header =TRUE, sep="", fill=TRUE)
ensembl.genes = input_data[[1]]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
    filters="ensembl_gene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id"),
    values=ensembl.genes,
    mart=mart)
genes

#save the data without N/A
colSums(is.na(genes)) 
##rowSums(is.na(genes))
genes_dropNA = genes[rowSums(is.na(genes))==0,]
save(genes_dropNA, file = "genes_dropNA.rda")
rm(genes_dropNA)

#merge with inputFile and convered gene id
#change column name
colnames(genes) = c("Ensembl", "entrez_id")
output = merge(x=input_data,
      y=genes,
      by='Ensembl',
      all=TRUE)
#save file
save(output, file = "conversion_to_entrezID.rda")
write.table(output, paste0(getwd(), "/conversion_to_entrezID.txt"), quote=FALSE, append=FALSE, row.names = FALSE)
