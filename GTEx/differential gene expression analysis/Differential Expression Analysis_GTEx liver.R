library(DESeq2)
library(tidyverse)
library("ggplot2")
library("apeglm") #Library for Log fold change shrinkage
library("ashr") #Library for Log fold change shrinkage
library("pheatmap")
library("ggplot2")

load(GTEx liver)
load(GTEx liver subject sex)
table(GTEx_sample_id$sex)

summary for GTEx data
female  male    NA  total 
2026     3410       226

#Reference for analysis method
http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#to prepare DESeq2 input file : countData, condition
countData = GTEx_liver
ncol(GTEx_liver)

condition = GTEx_sample_id$Sex 
length(condition)

##check
ncol(GTEx_liver) == length(condition)

###Differential Expression analysis (DESeq2)###
#To prepare
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = DataFrame(condition),
                              design = ~ condition) #generate the deseq data set

#filter
dds <- dds[ rowSums(counts(dds)) > 1, ] #remove genes with zero counts

#reference level
dds$condition #reference : female

#DEA analysis
dds <- DESeq(dds)
rm(countData, condition)

#DEA results
head(sizeFactors(dds))
res = results(dds)
summary(res)

#Log fold change shrinkage for visualization and ranking (3 type)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_male_vs_female", type="apeglm")
resLFC

#save the results
save(dds, res, resLFC, file='dds.rda')


#p-value and adjusted p-value
resOrdered <- res[order(res$pvalue),] #order our results table by the smallest p value
sum(res$padj < 0.1, na.rm=TRUE) #How many adjusted p-values were less than 0.1?

#if I want to change alpha value into 0.05,
# res05 <- results(dds, alpha=0.05)
# summary(res05)
# sum(res05$padj < 0.05, na.rm=TRUE)


#vsd <- vst(dds, blind = FALSE) #normalization considering tissue #I don't know 

##MA-plot
plotMA(res, main='MA-plot')
plotMA(res, ylim=c(-5,5), main='MA-plot(male vs female)')
plotMA(resLFC, ylim=c(-5,5), main='MA-plot (male vs female, apeglm)')

#Dispersion Plot
plotDispEsts(dds)

#More information on results columns
mcols(res)$description

#To export results to CSV files
write.csv(as.data.frame(resOrdered), 
          file="GTEx_liver_DEA-sex.csv")
resSig <- subset(resOrdered, padj < 0.1) #Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.
resSig

save(resOrdered, resSig, file='continue.rda')

#volcano plot
cut_lfc <- 1
cut_pvalue <- 0.01

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(resLFC)

## Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot (male_vs_female, apeglm)", col='grey', cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

rm(cut_lfc,cut_pvalue,topT)




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
save(vsd, sampleDists, file="vsd.rda")

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")

rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- GTEx_sample_id$sample_id
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#annotation
color_code=data.frame(sample = rownames(sampleDistMatrix))
rownames(color_code) = colnames(sampleDistMatrix)
#heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         cutree_rows = 2,
         cutree_cols = 2,
         col=colors,
         annotation_col = color_code,
         show_colnames = F,
         show_rownames = F,
         main = "Sample-to-sample distances_male vs female")


# Principal component plot of the samples
plotPCA(vsd, intgroup="condition")

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    ggtitle("PCA plot_male vs female")+
    coord_fixed()

rm(sampleDists, colors, percentVar, sampleDistMatrix, color_code, pcaData, vsd)
rm(resSig, resOrdered)

#extract resLFC list
dataframe_resLFC = data.frame(resLFC)
dataframe_resLFC_dropNA = drop_na(dataframe_resLFC)
male = dataframe_resLFC_dropNA[dataframe_resLFC_dropNA$log2FoldChange>1 & dataframe_resLFC_dropNA$padj>0.01,]
write.csv(male, 
          file="GTEx_liver_male.csv")
rm(dataframe_resLFC, dataframe_resLFC_dropNA, male)
