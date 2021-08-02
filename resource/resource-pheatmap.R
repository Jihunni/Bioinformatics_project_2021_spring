###pheatmap###
###draw Heatmap###

# Requirement1 : to create colors code for each group
##(ATP consuming group)
newCols_consuming <- colorRampPalette(grDevices::rainbow(length(unique(color_consuming$SUBSYSTEM)))) #color
annoCol_consuming <- newCols_consuming(length(unique(color_consuming$SUBSYSTEM)))
names(annoCol_consuming) <- unique(color_consuming$SUBSYSTEM) #name


##(ATP producing group)
rownames(color_producing) = color_producing$RXNID
newCols_producing <- colorRampPalette(grDevices::rainbow(length(unique(color_producing$SUBSYSTEM))))
annoCol_producing <- newCols_producing(length(unique(color_producing$SUBSYSTEM)))
names(annoCol_producing) <- unique(color_producing$SUBSYSTEM)

##(merge into one list)
annoCol <- list(subsystem_producing = annoCol_producing, subsystem_consuming = annoCol_consuming)


#Requirement2 : give annotation to give each data (rowname should match to name in correlation matix)
color_code_subsystem_consuming = data.frame(subsystem_consuming = color_consuming[,"subsystem_consuming"])
rownames(color_code_subsystem_consuming) = rownames(color_consuming)

color_code_subsystem_producing = data.frame(subsystem_producing = color_producing[,"subsystem_producing"])
rownames(color_code_subsystem_producing) = rownames(color_producing)




#draw heatmap
#for color annotation
colnames(color_consuming)[3] = 'subsystem_consuming'
colnames(color_producing)[3] = 'subsystem_producing'

pheatmap(mat = correlation_matrix[1:ncol(df1), as.numeric(ncol(df1)+1):ncol(correlation_matrix)],
         #color = colorRampPalette(c('#2471A3','white','#C0392B'))(50),
         dispaly_numbers = FALSE,
         border_color = 'white',
         cutree_rows = 2,
         cutree_cols = 2,
         annotation_colors = annoCol,
         annotation_row = color_code_subsystem_consuming, #it works
         #    annotation = color_code_subsystem_consuming, #Not works
         row_annotation_legend = TRUE, #It works
         
         annotation_col = color_code_subsystem_producing, #It works
         show_rownames = F,
         show_colnames = T,
         #    cellwidth = 1000,
         #    cellheight = 1000,
         main = "main title"
)


###ERROR###
Error in `[.data.frame`(color_consuming, , "subsystem_producing") : 
    undefined columns selected

-> You need to define the same row names for annotation and color_code_subsystem_consuming