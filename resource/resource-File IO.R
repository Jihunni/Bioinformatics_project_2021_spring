File I/O
#Prepare input file for DAVID
##foreground gene
load(dataframe_resLFC_dropNA)
foreground_geneList = substr(rownames(dataframe_resLFC_dropNA), start=1, stop=15)
write(foreground_geneList, 'David_propensity_foreground_geneList.txt')

##background gene
load(dataframe_resLFC)
background_geneList = substr(rownames(dataframe_resLFC), start=1, stop=15)
background_geneList = background_geneList[!background_geneList %in% foreground_geneList] 
write(background_geneList, 'David_propensity_background_geneList.txt')

#Input file for GSEA
rm(description,geneID_ensemble )
input_GSEA_normalized_count = counts(dds, normalized =TRUE)
rownames(input_GSEA_normalized_count) = substr(rownames(input_GSEA_normalized_count), start=1, stop=15)
description = data.frame(DESCRIPTION = rep(NA, nrow(input_GSEA_normalized_count)))
geneID_ensemble = data.frame(NAME = rownames(input_GSEA_normalized_count))
input_GSEA_normalized_count = cbind(geneID_ensemble, description, input_GSEA_normalized_count)
write.table(input_GSEA_normalized_count, "normalized count_GSEA input file_ensembleID.txt",quote=F,sep="\t", row.names = F)



##prepare the cls file for GSEA
library(plyr)
#cls = revalue(dds$condition, c("Female" = 0, "Male" = 1)) #change name of factor
cls = as.character(dds$condition)
cls = gsub("Female","1", cls)
cls = gsub("Male","2", cls)
cat("416 2 1", "\n", #num_sample ; num_class ; always 1
    "# Female Male ", "\n", #The order of the labels on the third line determines the association of class names and class labels
    cls, file = "cls.cls", append = FALSE)
#You need to remove space and check the seqeunce in GSEA program