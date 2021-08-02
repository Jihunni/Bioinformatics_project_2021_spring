library(readr)
library(stringr)
library(tidyverse)

#correlation
library("Hmisc") #correlation matrix with p-value
library(corrplot) #visualizing correlation matrix
library(pheatmap)
library(grid)
library(RColorBrewer) #color annotation

library(beepr) #beep() : alarm function

# #load pan-cancer data (small)
# tcga_tpm = read_tsv(file="../data/tcga_RSEM_gene_tpm",n_max=30)
# save(tcga_tpm, file = "../data/tcga_RSEM_gene_tpm_Only30.rda")

#load pan-cancer data (all)
tcga_tpm = read_tsv(file="../data/tcga_RSEM_gene_tpm")
save(tcga_tpm, file = "../data/tcga_RSEM_gene_tpm.rda")

<the numer of TCGA case per cancer type>
#load tissueSourceSite code data
tissueSourceSite = read_tsv(file="../tcga_code_tables/tissueSourceSite.tsv")
colname = colnames(tcga_tpm)
head(colname) #we need to consider the first column

#extract TSS code
TSS = substr(colname, start = 6, stop = 7) #tcga tpm data
head(TSS)

#HCC data.frame
TSS_code_HCC_boolian = str_detect(tissueSourceSite$`Study Name`, 'Liver hepatocellular carcinoma')
TSS_cide_HCC = tissueSourceSite$`TSS Code`[TSS_code_HCC_boolian]
Index_HCC = c(1, which(TSS %in% TSS_cide_HCC)) #sample column sholud be included
tcga_tpm_HCC = tcga_tpm[Index_HCC]
save(tcga_tpm_HCC, TSS_cide_HCC, Index_HCC, file = 'HCC_dataFrame.rda')
rm(tissueSourceSite, colname, Index_HCC, TSS, TSS_cide_HCC, TSS_code_HCC_boolian)


#sub-dataframe which contains ATP-consuming group and ATP-producing group
tcga_gene_list = substr(tcga_tpm_HCC$sample, start = 1, stop = 15) #ensemble ID has unqiue 11 digital
#gene list from HME data
load('../../GEM data/geneList_consuming.rda')
load('../../GEM data/geneList_producing.rda')

HCC_consuming = tcga_gene_list %in% geneList_consuming    
HCC_producing = tcga_gene_list %in% geneList_producing    
unique(HCC_consuming)
unique(HCC_producing)
HCC_consuming = tcga_tpm_HCC[c(HCC_consuming),]
HCC_producing = tcga_tpm_HCC[c(HCC_producing),]

save(HCC_consuming, HCC_producing, file='../data/HCC_consuming_and_producing_dataFrame')
rm(geneList_consuming, geneList_producing, mean_HCC, tcga_gene_list)


# #HCC gene expression median value
# plot(x = c(1:100), y = tcga_tpm_HCC[10,1:100], ylim=c(-15,15)) #to observe the overall trend
# mean_HCC_consuming = apply(HCC_consuming[,2:length(colnames(HCC_consuming))], 1, median)
# mean_HCC_producing = apply(HCC_producing[,2:length(colnames(HCC_producing))], 1, median)
# #tcga_tpm_HCC_mean = mutate(tcga_tpm_HCC, apply(tcga_tpm_HCC[,2:length(colnames(tcga_tpm_HCC))], 1,median))


# mean_HCC_consuming = data.frame(HCC_consuming$sample, mean_HCC_consuming)
# colnames(mean_HCC_consuming) = c('gene', 'tpm')
# mean_HCC_producing = data.frame(HCC_producing$sample, mean_HCC_producing)
# colnames(mean_HCC_producing) = c('gene', 'tpm')
# correlation_HCC = cor(mean_HCC_consuming$tpm, mean_HCC_producing$tpm)
# save(mean_HCC_consuming, mean_HCC_producing, file='HCC_mean_consuming_and_producing_dataFrame.rda')

<correlation analysis in gene level>
save(HCC_consuming, HCC_producing, tcga_tpm_HCC, file='continue.rda')


#Correlation matrix
##consuming group
row.names(HCC_consuming) = HCC_consuming$sample #change raw name into geneID
HCC_consuming_t = t(HCC_consuming) #transpose the dataframe (output: matrix)
HCC_consuming_t = HCC_consuming_t[2:nrow(HCC_consuming_t), ] #remove the first raw containing gene ID

##producing group
row.names(HCC_producing) = HCC_producing$sample #change raw name into geneID
HCC_producing_t = t(HCC_producing) #transpose the dataframe (output: matrix)
HCC_producing_t = HCC_producing_t[2:nrow(HCC_producing_t), ] #remove the first raw containing gene ID

##merge
unique(rownames(HCC_consuming_t) == rownames(HCC_producing_t))  #check whether they have identical raw
merged_matrix = merge(HCC_consuming_t, HCC_producing_t, by = "row.names", all = TRUE) #dataframe??? #merge two matrix
rownames(merged_matrix) = merged_matrix$Row.names
merged_matrix = merged_matrix[,2:ncol(merged_matrix)] #remove rowname column

is.na(merged_matrix) #check the data has N/A
which(is.na(merged_matrix))
merged_matrix_numberic = data.frame(lapply(merged_matrix, as.numeric)) #convert char dataframe into numeric dataframe.
rownames(merged_matrix_numberic) = rownames(merged_matrix)
merged_matrix = merged_matrix_numberic
rm(merged_matrix_numberic)

correlation_matrix = cor(merged_matrix)
save(HCC_consuming_t, HCC_producing_t, merged_matrix, correlation_matrix, file='correlation_matrix.rda')
load('correlation_matrix.rda')
#correlation analysis
Rcorrelation_matrix = rcorr(as.matrix(merged_matrix))
HCC_coeff = Rcorrelation_matrix$r
HCC_p = Rcorrelation_matrix$P

#correlogram
corrplot(correlation_matrix) #A default correlation matrix plot (called a Correlogram)

#Heatmap
palette = colorRampPalette(c("green", "white", "red")) (20) #Heatmap color
heatmap(x = correlation_matrix,
        col = palette, 
        symm = TRUE,
        RowSideColors = c(
            rep("red", ncol(HCC_consuming_t)), 
            rep("blue", ncol(HCC_producing_t))),
#        ColSideColors = c(
#            rep("red", ncol(HCC_consuming_t)), 
#            rep("blue", ncol(HCC_producing_t))),
        main = 'HCC correlation'
)


# ex = heatmap(x = correlation_matrix[1:10,1:10],
#         col = palette, 
#         symm = TRUE,
#         #RowSideColors = c(rep("red", 5), rep("blue", 5)),
#         main = 'HCC correlation'
# )


#color code
gene_color_code = data.frame(gene= rep(c("consuming", "producing"), c(ncol(HCC_consuming_t), ncol(HCC_producing_t))))
row.names(gene_color_code) = colnames(merged_matrix)

#pheatmap
# small_pheatmap_HCC = pheatmap(mat = correlation_matrix[1:10,1:10],
#                         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
#                         dispaly_numbers = FALSE,
#                         border_color = 'white',
#                         show_rownames = T,
#                         show_colnames = F,
#                         annotation_row = gene_color_code,
#                         annotation_col = gene_color_code
# 
# )

pheatmap_HCC = pheatmap(mat = correlation_matrix,
                        color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
                        dispaly_numbers = FALSE,
                        border_color = 'white',
                        show_rownames = T,
                        show_colnames = F,
                        annotation_row = gene_color_code,
                        annotation_col = gene_color_code 
#                        labels_row 
                        )


#pheatmap_HCC_consuming_vs_producing
#To prepare the color code for consuming group
gene_to_subsystem_consuming = colnames(HCC_consuming_t)
gene_to_subsystem_consuming = data.frame(gene = substr(gene_to_subsystem_consuming, start = 1, stop = 15))
unique(duplicated(gene_to_subsystem_consuming)) #no duplication
gene_to_subsystem_consuming = gene_to_subsystem_consuming %>% 
    left_join(y= gene_to_HMRdata_df, by = ('gene' = 'gene')) %>%
    select(gene, SUBSYSTEM)
gene_to_subsystem_consuming = distinct(gene_to_subsystem_consuming)
    #eventhough this function remove duplicated one, it would not change sequence of gene list. Look at the example of left_join

##issue : dupliated element owing to multiple matching to subsystem
###multiple matching --> change into 'multiple'
1. find out duplicated
dupliated_consuming_list = gene_to_subsystem_consuming[duplicated(gene_to_subsystem_consuming$gene),]

2. remove duplicated list and make it one
gene_to_subsystem_consuming = gene_to_subsystem_consuming[!(duplicated(gene_to_subsystem_consuming$gene)),]

3. changed the duplicated list into 'multiple'
gene_to_subsystem_consuming[gene_to_subsystem_consuming$gene %in% dupliated_consuming_list$gene,]$SUBSYSTEM = 'multiple'




gene_to_subsystem_producing = colnames(HCC_producing_t)
gene_to_subsystem_producing = data.frame(gene = substr(gene_to_subsystem_producing, start = 1, stop = 15))
gene_to_subsystem_producing = gene_to_subsystem_producing %>% 
    left_join(y= gene_to_HMRdata_df, by = ('gene' = 'gene')) %>%
    select(gene, SUBSYSTEM)
gene_to_subsystem_producing = distinct(gene_to_subsystem_producing)

##same issue : dupliated element owing to multiple matching to subsystem
###multiple matching --> change into 'multiple'
1. find out duplicated
dupliated_producing_list = gene_to_subsystem_producing[duplicated(gene_to_subsystem_producing$gene),]

2. remove duplicated list and make it one
gene_to_subsystem_producing = gene_to_subsystem_producing[!(duplicated(gene_to_subsystem_producing$gene)),]

3. changed the duplicated list into 'multiple'
gene_to_subsystem_producing[gene_to_subsystem_producing$gene %in% dupliated_producing_list$gene,]$SUBSYSTEM = 'multiple'

save(HCC_consuming_t, HCC_producing_t, gene_to_HMRdata_df, gene_to_subsystem_consuming, gene_to_subsystem_producing, dupliated_consuming_list, dupliated_producing_list, file = 'HCC_color_code.rda')

#color code for subsystem
color_code_subsystem_consuming = data.frame(subsystem = gene_to_subsystem_consuming$SUBSYSTEM)
color_code_subsystem_producing = data.frame(subsystem = gene_to_subsystem_producing$SUBSYSTEM)

#final check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
unique(substr(colnames(HCC_consuming_t), start = 1, stop = 15) == gene_to_subsystem_consuming$gene)
unique(substr(colnames(HCC_producing_t), start = 1, stop = 15) == gene_to_subsystem_producing$gene)
    #GOOOOOOOOOOOOOOD
rownames(color_code_subsystem_consuming) = colnames(HCC_consuming_t)
rownames(color_code_subsystem_producing) = colnames(HCC_producing_t)

ex_colorcode_col = data.frame(subsystem = rep(1, ncol(HCC_consuming_t)))
ex_colorcode_row = data.frame(subsystem = rep(2, ncol(correlation_matrix) - ncol(HCC_consuming_t)))

#to escape from collision
colnames(color_code_subsystem_consuming) = 'subsystem_consuming'    
colnames(color_code_subsystem_producing) = 'subsystem_producing'

# to create colors for each group 
##(ATP consuming group)
newCols_consuming <- colorRampPalette(grDevices::rainbow(length(unique(color_code_subsystem_consuming$subsystem_consuming))))
annoCol_consuming <- newCols_consuming(length(unique(color_code_subsystem_consuming$subsystem_consuming)))
names(annoCol_consuming) <- unique(color_code_subsystem_consuming$subsystem_consuming)
##(ATP producing group)
newCols_producing <- colorRampPalette(grDevices::rainbow(length(unique(color_code_subsystem_producing$subsystem_producing))))
annoCol_producing <- newCols_producing(length(unique(color_code_subsystem_producing$subsystem_producing)))
names(annoCol_producing) <- unique(color_code_subsystem_producing$subsystem_producing)

annoCol <- list(subsystem_consuming = annoCol_consuming, subsystem_producing = annoCol_producing)

#draw heatmap
pheatmap(mat = correlation_matrix[1:ncol(HCC_consuming_t), as.numeric(ncol(HCC_consuming_t)+1):ncol(correlation_matrix)],
        color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
        dispaly_numbers = FALSE,
        border_color = 'white',
        cutree_rows = 2,
        cutree_cols = 2,
        annotation_colors = annoCol,
        annotation_row = color_code_subsystem_consuming, #it works
    #    annotation = color_code_subsystem_consuming, #Not works
    #   row_annotation_legend = FALSE, #It works
        
        annotation_col = color_code_subsystem_producing, #It works
        show_rownames = F,
        show_colnames = F,
    #    cellwidth = 1000,
    #    cellheight = 1000,
        main = "HCC (x: ATP producing_192, y: ATP consuming_1084)"
        )



 
save(color_code_subsystem_consuming, color_code_subsystem_producing, annoCol_consuming, annoCol_producing, correlation_matrix, HCC_consuming_t, HCC_producing_t, file='continue.rda')
#this figure sholud be exported in big size due to labeling

# alarm()
# hclust_HCC = hclust(dist())

save(correlation_matrix, Rcorrelation_matrix, HCC_coeff, HCC_p, file='HCC_correlation.rda')

#plot(c(1,2), c(1,2)) #check RStuido deos work right now.

<Further analysis for oxidative phospohrylation>
producing_oxidative_phosphorylation_df = gene_to_subsystem_producing %>%
    left_join(y=gene_to_HMRdata_df, by = c('gene' ='gene')) %>%
    select(gene, RXNID, EQUATION, SUBSYSTEM.x)
producing_oxidative_phosphorylation_table = table(producing_oxidative_phosphorylation_df$RXNID)
plot(x=factor(c('HMR_4358', 'HMR_4421', 'HMR_4572', 'HMR_4675', 'HMR_5295', 'HMR_5429', 'HMR_6916', 'HMR_7629', 'HMR_7799', 'HMR_8474', 'HMR_8892', 'HMR_9570')), y=c(4,2,9,9,5,5,178,9,10,5,3,1), ylim=c(0,200))

###filtering on Heatmap###
#filter on ATP consuming group :  #gene >10 per subsystem
load(filter_subsystem_consuming)
filter_gene_to_subsystem_consuming = gene_to_subsystem_consuming[gene_to_subsystem_consuming$SUBSYSTEM %in% filter_subsystem_consuming$SUBSYSTEM,]

filter_gene_boolian = c((substr(colnames(merged_matrix), start = 1, stop = 15) %in% filter_gene_to_subsystem_consuming$gene)[1:1084], rep(TRUE, nrow(gene_to_subsystem_producing)))


filter_merged_matrix =  merged_matrix[,filter_gene_boolian]
filter_correlation_matrix = correlation_matrix[filter_gene_boolian,filter_gene_boolian]

#annotation#
color_code_subsystem_consuming = data.frame(gene = gene_to_subsystem_consuming$gene, subsystem_consuming = gene_to_subsystem_consuming$SUBSYSTEM)
color_code_subsystem_consuming = color_code_subsystem_consuming[color_code_subsystem_consuming$subsystem %in% filter_subsystem_consuming$SUBSYSTEM,] #filtering
rownames(color_code_subsystem_consuming) = color_code_subsystem_consuming$gene
color_code_subsystem_consuming = select(color_code_subsystem_consuming, subsystem_consuming)

color_code_subsystem_producing = data.frame(subsystem_producing = gene_to_subsystem_producing$SUBSYSTEM)
rownames(color_code_subsystem_producing) = gene_to_subsystem_producing$gene

# to create colors for each group 
##(ATP consuming group)
newCols_consuming <- colorRampPalette(grDevices::rainbow(length(unique(color_code_subsystem_consuming$subsystem_consuming))))
annoCol_consuming <- newCols_consuming(length(unique(color_code_subsystem_consuming$subsystem_consuming)))
names(annoCol_consuming) <- unique(color_code_subsystem_consuming$subsystem_consuming)

##(ATP producing group)
newCols_producing <- colorRampPalette(grDevices::rainbow(length(unique(color_code_subsystem_producing$subsystem_producing))))
annoCol_producing <- newCols_producing(length(unique(color_code_subsystem_producing$subsystem_producing)))
names(annoCol_producing) <- unique(color_code_subsystem_producing$subsystem_producing)

annoCol <- list(subsystem_consuming = annoCol_consuming, subsystem_producing = annoCol_producing)


#final check for annotation
unique(substr(rownames(filter_correlation_matrix)[1:896], start = 1, stop = 15) == color_code_subsystem_consuming$gene) #It sholud be TRUE
unique(substr(rownames(filter_correlation_matrix)[897:1088], start = 1, stop = 15) == rownames(color_code_subsystem_producing)) #It sholud be TRUE

#to fit the dataframe format in heatmap requirement
rownames(color_code_subsystem_consuming) = rownames(filter_correlation_matrix)[1:896]
rownames(color_code_subsystem_producing) = rownames(filter_correlation_matrix)[897:1088]
   
#draw heatmap
pheatmap(mat = filter_correlation_matrix[1:896, 897:1088],
         color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         dispaly_numbers = FALSE,
         border_color = 'white',
         cutree_rows = 2,
         cutree_cols = 2,
         annotation_colors = annoCol,
         annotation_row = color_code_subsystem_consuming, #it works
         annotation_col = color_code_subsystem_producing, #It works
         show_rownames = F,
         show_colnames = F,
         #    cellwidth = 1000,
         #    cellheight = 1000,
         main = "HCC with filter (x: ATP producing_192, y: ATP consuming_896)"
)

save(color_code_subsystem_consuming, color_code_subsystem_producing, filter_correlation_matrix, filter_gene_boolian, filter_gene_to_subsystem_consuming, filter_merged_matrix, filter_subsystem_consuming, gene_to_subsystem_consuming, gene_to_subsystem_producing, file='filter_HCC_heatmap.rda')
