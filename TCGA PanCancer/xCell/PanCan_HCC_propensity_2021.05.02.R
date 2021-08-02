library(tidyverse)
library(dplyr) 
library(ggplot2)
library("ggpubr") #ggplot with p-value

#load files
load(tcga_RSEM_gene_TPM)
load(propensity_normalized)

HCC_tpm = tcga_tpm[propensity_normalized$case]
rownames(HCC_tpm) = rownames(tcga_tpm)
save(HCC_tpm, file='PanCan_HCC_tpm.rda')
example_input = read_tsv('example_input.txt')
gene_list = read_tsv('gene_list.txt', col_names = FALSE)

rownames(HCC_tpm) = gene_list$X1

#annotation, biomart
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- gene_list$X1
ensLookup <- gsub("\\.[0-9]*$", "", ens)
ensLookup

annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)

annotLookup <- data.frame(
    ens[match(annotLookup$ensembl_gene_id, ensLookup)],
    annotLookup)

colnames(annotLookup) <- c(
    "original_id",
    c("ensembl_gene_id", "gene_biotype", "external_gene_name"))

annotLookup

colnames(annotLookup)[1] ='gene_id'
gene_list = left_join(gene_list,annotLookup, by=c('X1'='gene_id'))
rm(ens, ensLookup,mart)


#merge with count file and annotation file
output = HCC_tpm[!is.na(gene_list$external_gene_name),]
output_rownames = gene_list$external_gene_name[!is.na(gene_list$external_gene_name)]
output = cbind(output_rownames, output)

write.table(output, file = "xCell_input_PanCan_HCC_2021.05.02.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE) 


#propensity algorithm
colnames(output) == propensity_normalized$case
weight_vector = data.frame(case=colnames(output)[-1])
vector = propensity_normalized[,c(1,11)]
weight_vector = left_join(weight_vector, vector, by=c("case"="case"))
weight_vector = as.vector(weight_vector)

output_propensity = data.frame(mapply('*',output[,-1],weight_vector))
output_propensity = cbind(output_rownames, output_propensity)

write.table(output_propensity, file = "xCell_input_PanCan_HCC_propensity_2021.05.02.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE) 

save(gene_list, output, output_propensity, propensity_normalized, file='2021.05.02.rda')
rm(tcga_tpm, vector, weight_vector, example_input, HCC_tpm)
rm(annotLookup, output_rownames)


#read the output from xCell
xCell_result = read_tsv('xCell_output_PanCan_HCC_propensity_2021.05.02_final result.txt')
xCell_pval = read_tsv('xCell_output_PanCan_HCC_propensity_2021.05.02_pvals.txt')

colnames(xCell_result) = gsub("[.]","-", colnames(xCell_result)) #change colname

save(xCell_result, xCell_pval, file='xcell_output_2021.05.03.rda') #save output file


### plot ###

#grouping in propensity normalized file
#Tumor and TumorFree in Male
male_tumor = propensity_normalized[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'With Tumor',]$case

male_tumorFree = propensity_normalized[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'Tumor Free',]$case

# table(colnames(HCC_protein) %in% male_tumor) #24
# table(colnames(HCC_protein) %in% male_tumorFree) #82
# HCC_protein[,!colnames(HCC_protein) %in% male_tumor & !colnames(HCC_protein) %in% male_tumorFree]

#Tumor w/o hepatitis, Hepatitis B, Hepatitis C in male
male_tumor_no = propensity_normalized[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'With Tumor' & propensity_normalized$hepatitis=='no',]$case
male_tumor_B = propensity_normalized[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'With Tumor' & propensity_normalized$hepatitis=='B',]$case
male_tumor_C = propensity_normalized[propensity_normalized$sex == 'Male' & propensity_normalized$status == 'With Tumor' & propensity_normalized$hepatitis=='C',]$case

# table(colnames(HCC_protein) %in% male_tumor_no) #15
# table(colnames(HCC_protein) %in% male_tumor_B) #3
# table(colnames(HCC_protein) %in% male_tumor_C) #5

#Tumor and TumorFree in Male
female_tumor = propensity_normalized[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'With Tumor',]$case

female_tumorFree = propensity_normalized[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'Tumor Free',]$case

# table(colnames(HCC_protein) %in% female_tumor) #25
# table(colnames(HCC_protein) %in% female_tumorFree) #39
# HCC_protein[,!colnames(HCC_protein) %in% female_tumor & !colnames(HCC_protein) %in% female_tumorFree]

#Tumor w/o hepatitis, Hepatitis B, Hepatitis C in female
female_tumor_no = propensity_normalized[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'With Tumor' & propensity_normalized$hepatitis=='no',]$case
female_tumor_B = propensity_normalized[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'With Tumor' & propensity_normalized$hepatitis=='B',]$case
female_tumor_C = propensity_normalized[propensity_normalized$sex == 'Female' & propensity_normalized$status == 'With Tumor' & propensity_normalized$hepatitis=='C',]$case

# table(colnames(HCC_protein) %in% female_tumor_no) #23
# table(colnames(HCC_protein) %in% female_tumor_B) #1
# table(colnames(HCC_protein) %in% female_tumor_C) #1

# Violin plots with box plots inside - for loop
##detach("package:dplyr", unload=TRUE)
for(i in xCell_result$type){
#    i = 'Astrocytes'
    print(i)
    gene = i
    index = which(xCell_result$type == gene)
    #all group comparison between male and female
    # gene_male_tumorFree %>% 
    #     select(colnames(xCell_result)[colnames(xCell_result) %in% male_tumorFree]) %>%
    #     filter()
    #     as.numeric(xCell_result[gene, colnames(xCell_result) %in% male_tumorFree])
    gene_male_tumorFree = as.numeric(xCell_result[index, colnames(xCell_result) %in% male_tumorFree])
    gene_male_tumor = as.numeric(xCell_result[index, colnames(xCell_result) %in% male_tumor])
    gene_male_tumor_no = as.numeric(xCell_result[index, colnames(xCell_result) %in% male_tumor_no])
    gene_male_tumor_B = as.numeric(xCell_result[index, colnames(xCell_result) %in% male_tumor_B])
    gene_male_tumor_C = as.numeric(xCell_result[index, colnames(xCell_result) %in% male_tumor_C])
    
    
    gene_female_tumorFree = as.numeric(xCell_result[index, colnames(xCell_result) %in% female_tumorFree])
    gene_female_tumor = as.numeric(xCell_result[index, colnames(xCell_result) %in% female_tumor])
    
    ##mean(gene_male)
    ##mean(gene_female)
    ##t.test(gene_male, gene_female)
    
    #input data for plot
    gene_data = data.frame(value=c(gene_male_tumorFree, gene_male_tumor, gene_male_tumor_no, gene_male_tumor_B, gene_male_tumor_C, gene_female_tumorFree, gene_female_tumor), 
                           type=c(rep('male_tumorFree', length(gene_male_tumorFree)), rep('male_tumor', length(gene_male_tumor)), rep('male_tumor_no', length(gene_male_tumor_no)), rep('male_tumor_B', length(gene_male_tumor_B)), rep('male_tumor_C', length(gene_male_tumor_C)), rep('female_tumorFree', length(gene_female_tumorFree)), rep('female_tumor', length(gene_female_tumor))))
    
    gene_data = mutate(gene_data, model=factor(type, levels=c("female_tumorFree", 'female_tumor', 'male_tumorFree', 'male_tumor', 'male_tumor_no', 'male_tumor_B', 'male_tumor_C')))
    
    #violin plot with box plot
    comparison=list(c("female_tumorFree", "female_tumor"), c("male_tumorFree", "male_tumor"), c("female_tumor", "male_tumor"),c("male_tumor_no", "male_tumor_B"), c("male_tumor_no", "male_tumor_C"))
    
    ggviolin(gene_data, x = "model", y = "value", fill = "model",
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"))+
        stat_compare_means(aes(type), comparisons = comparison) + # Add significance levels
    #    +stat_compare_means(label.y = 3)        # Add global the p-value 
        ggtitle(gene)  + ylab("abundance")
    
    ggsave(paste0("./figure/",as.character(gene),".png"), width=25, height = 15, units='cm', limitsize = FALSE)
    rm(comparison, gene, gene_data)
    rm(gene_male_tumor, gene_male_tumor_B, gene_male_tumor_C, gene_male_tumor_no, gene_male_tumorFree, gene_female_tumor, gene_female_tumorFree)
} #before running this code, check the folder file
