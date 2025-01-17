package

#library
library(readr)
library(ggplot2)
library(tidyverse)
library(readxl)
library("Hmisc") #correlation matrix with p-value
library(corrplot) #visualizing correlation matrix
library(beepr) #beep() : alarm function
library(DESeq2)
library("apeglm") #shrink log fold change
library(ggraph)
library(igraph)
library("vsn")
library(piano) #report metabolite analysis
#install package
require()

packageVersion("TCGAbiolinks")

install.packages('ggraph')

if (!"BiocManager" %in% rownames(installed.packages()))
    install.packages("dendextend")
BiocManager::install("piano")

sessionInfo()

##Error and trials
ERROR: failed to lock directory ‘C:\R\R-4.0.3\library’ for modifying
Try removing ‘C:\R\R-4.0.3\library/00LOCK
1. install.packages("설치하고 싶은 패키지 이름", dependencies = TRUE, INSTALL_opts = "--no-lock")
2. find the file:‘C:\Program Files\R\R-3.6.1\library/00LOCK’, delete the file"00LOCK" then it works again


###Bioinformatics###
#convert gene name to ensemble gene ID
library("AnnotationDbi")
library("org.Hs.eg.db")
options(connectionObserver = NULL)
gene_id = mapIds(org.Hs.eg.db,
                 keys=TCGA_LIHC_gene_count$sample, 
                 column="ENSEMBL",
                 keytype="SYMBOL",
                 multiVals="first")
gene_id = data.frame(gene_id)
nrow(gene_id)
length(which(gene_id$gene_id != 'NA')) #17311
TCGA_LIHC_gene_count = cbind(gene_id, TCGA_LIHC_gene_count) #add sample column at first


#annotation, biomart
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listMarts()
mart <- useDataset("hsapiens_gene_ensembl", mart)
listDatasets(mart)
filter = listFilters(mart)

ens <- c("ENSG00000100601.5", "ENSG00000178826.6",
         "ENSG00000243663.1", "ENSG00000138231.8")
ensLookup <- gsub("\\.[0-9]*$", "", ens)

annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_transcript_id", "ensembl_gene_id",
                 "gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)

annotLookup <- data.frame(
    ens[match(annotLookup$ensembl_gene_id, ensLookup)],
    annotLookup)

colnames(annotLookup) <- c(
    "original_id",
    c("ensembl_transcript_id", "ensembl_gene_id",
      "gene_biotype", "external_gene_name"))

#anntoation biomart version 2 (no digital below 0)
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listMarts()
mart <- useDataset("hsapiens_gene_ensembl", mart)
listDatasets(mart)
filter = listFilters(mart)

ensLookup <- gene_id
annotLookup <- getBM(
  mart=mart,
  attributes=c(,"ensembl_gene_id","gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensLookup,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensLookup[match(annotLookup$ensembl_gene_id, ensLookup)], annotLookup)

colnames(annotLookup) <- c(
  "original_id",
  c("ensembl_transcript_id", "ensembl_gene_id",
    "gene_biotype", "external_gene_name"))