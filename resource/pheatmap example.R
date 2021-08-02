#example of pheatmap

# Generate some data
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# Draw heatmaps
pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)
pheatmap(test, display_numbers = TRUE)
pheatmap(test, display_numbers = TRUE, number_format = "%.1e")
pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
                                                                            "1e-4", "1e-3", "1e-2", "1e-1", "1"))
pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
#pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")


# Generate column annotations
annotation = data.frame(Var1 = factor(1:10 %% 2 == 0,
                                      labels = c("Class1", "Class2")), Var2 = 1:10)
annotation$Var1 = factor(annotation$Var1, levels = c("Class1", "Class2", "Class3"))
rownames(annotation) = paste("Test", 1:10, sep = "")

pheatmap(test, annotation = annotation)
pheatmap(test, annotation = annotation, annotation_legend = FALSE)
pheatmap(test, annotation = annotation, annotation_legend = FALSE, drop_levels = FALSE)

# Specify colors
Var1 = c("navy", "darkgreen")
names(Var1) = c("Class1", "Class2")
Var2 = c("lightgreen", "navy")

ann_colors = list(Var1 = Var1, Var2 = Var2)

#Specify row annotations
row_ann <- data.frame(foo=gl(2,nrow(test)/2),`Bar`=relevel(gl(2,nrow(test)/2),"2"))
rownames(row_ann)<-rownames(test)
pheatmap(test, annotation = annotation, annotation_legend = FALSE, drop_levels = FALSE,row_annotation = row_ann)

#Using cytokine annotations
M<-matrix(rnorm(8*20),ncol=8)
row_annotation<-data.frame(A=gl(4,nrow(M)/4),B=gl(4,nrow(M)/4))
eg<-expand.grid(factor(c(0,1)),factor(c(0,1)),factor(c(0,1)))
colnames(eg)<-c("IFNg","TNFa","IL2")
rownames(eg)<-apply(eg,1,function(x)paste0(x,collapse=""))
rownames(M)<-1:nrow(M)
colnames(M)<-rownames(eg)
cytokine_annotation=eg
pheatmap(M,annotation=annotation,row_annotation=row_annotation,annotation_legend=TRUE,row_annotation_legend=TRUE,cluster_rows=FALSE,cytokine_annotation=cytokine_annotation,cluster_cols=FALSE)

# Specifying clustering from distance matrix
drows = dist(test, method = "minkowski")
dcols = dist(t(test), method = "minkowski")
pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)