#Convert gene name from TCGA data to Ensemble gene id

#check the library
require("mygene")

#install the library
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mygene")

# do the job (example)
library(mygene)

mygenes=c("TP53","AGTR1")
result = queryMany(mygenes, scopes="symbol", fields="ensembl.gene", species="human")
result[,c(1,4)]

# do the real job
input_data = read.table("C:/R/default_working_directory/Secretome/BRCA Project/output/dataDEGs.txt", header =TRUE, sep="")
gene_name = row.names(input_data)
result = queryMany(gene_name, scopes="symbol", fields="ensembl.gene", species="human")
result

write.table(result, paste0(getwd(), "/conversion.txt"), quote=FALSE, append=FALSE)


# count the number of N/A
input_data = read.table(paste0(getwd(), "/conversion.txt"), header =TRUE, sep="")
input_data = read.table(paste0(getwd(), "/conversion.txt"), header =FALSE, sep="", skip=1, skipNul = TRUE)
