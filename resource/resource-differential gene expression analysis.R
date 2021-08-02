library(readr)
library(tidyverse)
library(DESeq2)
library("ggplot2")
library("apeglm") #Library for Log fold change shrinkage
##library("ashr") #Library for Log fold change shrinkage
library("pheatmap")
library("ggplot2")

#Reference for analysis method
http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#merge GTEx and TCGA (HCC)
load(GTEx_liver)

#HCC without hepatitis
tumor_noHepatitis.propensity.data = tumor_noHepatitis.propensity.data[tumor_noHepatitis.propensity.data$weights != 0, ]
HCC_geneCount_in_Clinic = HCC_geneCount[tumor_noHepatitis.propensity.data$case]
rownames(HCC_geneCount_in_Clinic) = substr(rownames(HCC_geneCount), start=1, stop=15)

HCC_geneCount_in_Clinic = 2^HCC_geneCount_in_Clinic - 1 #do not forget!


table(tumor_noHepatitis.propensity.data$sex)
Female   Male   total
33     32       65

#summary
sample number   gene number
HCC_noH     65              60,498
GTEx_liver  226             56,200

#the number of matched gene
match = substr(rownames(HCC_geneCount_in_Clinic), start = 1, stop = 15) %in% rownames(GTEx_liver)
length(match[match == TRUE]) #the number of matched gene is 55475
rm(match)
#IT's good number to go further toward the next analysis


#merge two data into one data.frame
##data.frame(sample = rownames(GTEx_liver))
GTEx_liver = cbind(data.frame(sample = rownames(GTEx_liver)), GTEx_liver) #add sample column at first
HCC_geneCount_in_Clinic = cbind(data.frame(sample = rownames(HCC_geneCount_in_Clinic)), HCC_geneCount_in_Clinic)#add sample column at first

merged_df = inner_join(x = GTEx_liver, y = HCC_geneCount_in_Clinic, by = c('sample' = 'sample'))
save(merged_df, file='GTEx_and_HCCclinic_2021.05.15.rda')

#to prepare DESeq2 input file  : countData, condition
countData = merged_df[,-c(1)]
countData = data.frame(lapply(countData, as.integer)) #convert numeric datatype into integer type
rownames(countData) = merged_df$sample
rm(merged_df)

##summary
#summary
sample number   gene number
HCC_noH     65              60,498
GTEx_liver  226             56,200
merged_df   291             55,475      GTEx(1:226), TCGA(227,291)

condition = factor(c(rep('normal_liver',226), rep('HCC',65)))

save(countData, condition, file='DESeq2_input_2021.05.15.rda')
rm(merged_df, GTEx_liver, tcga_geneCount_HCC)

###Differential Expression analysis (DESeq2)###
#To prepare
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = DataFrame(condition),
                              design = ~ condition) #generate the deseq data set
#filter
dds <- dds[ rowSums(counts(dds)) > 1, ] #remove genes with zero counts

#Note on filter level
#dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
dds$condition <- relevel(dds$condition, ref = "normal_liver") #to specify the reference level
#dds$condition <- droplevels(dds$condition)

#DEA analysis
dds <- DESeq(dds)
rm(countData, condition)
#DEA results
head(sizeFactors(dds))
res = results(dds) #default alpha is 0.1
summary(res)

#Log fold change shrinkage for visualization and ranking (3 type)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_HCC_vs_normal_liver", type="apeglm")

save(dds, resLFC, file='dds_result_2021.05.15.rda')
resultsNames(dds)


#Exploring and exporting results
##MA-plot
plotMA(resLFC, ylim=c(-12,12), main='MA-plot (apeglm)')
plotMA(resLFC, ylim=c(-1,1), main='MA-plot (apeglm, small)')

idx <- identify(res$baseMean, res$log2FoldChange) #to interactively detect the row number of individual genes by clicking on the plot 
rownames(res)[idx] #recover the gene identifiers by saving the resulting indices

par(mfrow=c(1,4), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="None")
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
#mtext("MA-plot (small scale)", side = 3, line = 5, outer = TRUE)

#Plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition", main='Plot Counts')

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks=c(25,100,400)) +
    ggtitle('Plot Counts (smallest p-value gene)')

#volcano plot
cut_lfc <- 1
cut_pvalue <- 0.01

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC)

## Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot (HCC tumor without viral hepatitis vs normal liver tissue)", col='grey', cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

rm(cut_lfc,cut_pvalue,topT, resOrdered)


#Dispersion Plot
plotDispEsts(dds)

#More information on results columns
mcols(res)$description




### Data transformations and visualization ###
# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

# Effects of transformations on the variance
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
head(assay(vsd), 3)

meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")

rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- colnames(countData)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#annotation
color_code=data.frame(sample = rownames(sampleDistMatrix))
rownames(color_code) = colnames(sampleDistMatrix)
#heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         #  cutree_rows = 2,
         #  cutree_cols = 2,
         col=colors,
         annotation_col = color_code,
         show_colnames = F,
         show_rownames = F,
         main = "Sample-to-sample distances(HCC without hepatitis))")


# Principal component plot of the samples
##plotPCA(vsd, intgroup="condition")

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    ggtitle("PCA plot_PanCan(HCC without viral hepatitis)")+
    coord_fixed()



#convert gene name to ensemble gene symbol
df = data.frame(resLFC)
df = cbind(gene=rownames(df), df)
rownames(df) = 1:nrow(df)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

ens = ensLookup = df$gene
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

table(is.na(annotLookup)) #check NA
annotLookup = annotLookup[, c("ensembl_gene_id", "gene_biotype", "external_gene_name")]
df = left_join(df, annotLookup, by=c("gene" = "ensembl_gene_id"))
df = df[,c(1,7,8,1,2,3,4,5,6)]

rm(mart, annotLookup, ens, ensLookup)
save(df, file='dds_df_2021.05.15.rda')
write.table(df, file = "HCC_vs_GTEx_2021.05.15.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE)

df_dropNA = drop_na(df, padj, lfcSE)
write.table(df_dropNA, file = "HCC_vs_GTEx_dropNA_2021.05.15.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE)

table(df$pvalue < 0.05 & df$log2FoldChange > 1)


#GSEA
#Input file for GSEA
rm(description,geneID_ensemble )
input_GSEA_normalized_count = counts(dds, normalized =TRUE)
rownames(input_GSEA_normalized_count) = substr(rownames(input_GSEA_normalized_count), start=1, stop=15)
description = data.frame(DESCRIPTION = rep(NA, nrow(input_GSEA_normalized_count)))
geneID_ensemble = data.frame(NAME = rownames(input_GSEA_normalized_count))
input_GSEA_normalized_count = cbind(geneID_ensemble, description, input_GSEA_normalized_count)
write.table(input_GSEA_normalized_count, "GSEA input file_normalizedCount_ensembleID_2021.05.15.txt",quote=F,sep="\t", row.names = F)

rm(description,geneID_ensemble )



##prepare the cls file for GSEA
library(plyr)
#cls = revalue(dds$condition, c("Female" = 0, "Male" = 1)) #change name of factor
cls = as.character(dds$condition)

unique(dds$condition)
table(dds$condition)
unique(cls)
cat("291 2 1", "\n", #num_sample ; num_class ; always 1
    "# normal_liver HCC", "\n", #The order of the labels on the third line determines the association of class names and class labels
    cls, file = "GSEA_input_2021.05.15.cls", append = FALSE)

rm(cls, input_GSEA_normalized_count)


#visualize pathway
library(pathview)
data(gse16873.d)
data(demo.paths)
data(gene.idtype.list)
demo.paths$sel.paths
pathview_data = df$log2FoldChange
names(pathview_data) = rownames(df)
pathview(gene.data = pathview_data, gene.idtype ='ENSEMBL', pathway.id = "04630", species = "human", kegg.dir="C:/R/default_working_directory/Data/KEGG_pathway/", out.suffix = 'JAK-STAT signaling pathway')