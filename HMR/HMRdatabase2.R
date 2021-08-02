library(readxl)
library(tidyverse)

#load database
HMR = read_excel(path = './HMRdatabase2_00.xlsx', 
                 sheet = 'RXNS',# sheet name to read from
                 range = 'B1:S8187', # cell range to read from
                 col_names = TRUE,# TRUE to use the first row as column names
                 col_type = "guess", # guess the types of columns
                 na = "NA") # Character vector of strings to use for missing values
save(HMR, file = "HMR_DataFrame.rda")

#detect mathces
f = function(input_char){
    return(str_detect(input_char, "ATP"))
}

result = sapply(HMR$EQUATION, f)
    #HMR$EQUATION : vector
save(HMR, file = "HMR_ATP_detect.rda")

HMR_2 = select(HMR, c('RXNID', 'EQUATION', 'EC-NUMBER', 'GENE ASSOCIATION', 'COMPARTMENT', 'SUBSYSTEM'))
save(HMR_2, file = "HMR_2.rda")

HMR_3 = mutate(HMR_2, result)
save(HMR_3, file = "HMR_3.rda")

#Quesiton : multiple matches? 





#make a list which contains ATP
HMR_3_filter = filter(HMR_3, result == TRUE)
HMR_3_filter = HMR_3_filter[-c(7)]
save(HMR_3_filter, file = "HMR_3_filter.rda")
write.table(HMR_3_filter, file = "HNR_ATP_filter.txt", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#make a gene list which contains ATP
geneList = select(HMR_3_filter, 'GENE ASSOCIATION')
geneList = drop_na(geneList)

#geneList_2 = str_split_fixed(geneList$`GENE ASSOCIATION`, ';') 
geneList_2 = strsplit(geneList$`GENE ASSOCIATION`, ';')  #convert dataframe into list 
geneList_3 = unlist(geneList_2) # convert list into vector
#write.table(geneList_3, file = "geneList.txt", 
#            sep = "\t", col.names = FALSE, row.names = FALSE)
save(geneList_3, file = "geneList_vector.rda")
write(geneList_3, "geneList_3.txt")


<To separate Reaction Equation into three Group>
#Process
1.separate equation into two group. One is uni-directional, and the other is bi-directional
2. Among uni-directional group, separate left and right by '=>' using strsplit function
3. Get boolian result table (left group, right group)
4. get factor table. 1 stands for ATP-consuming. 2 stands fir ATP-producing


1. separate equation into two group
#separate dataframe into equilbirum group and non-equilibrum group
equilibrum_index = str_detect(HMR_3_filter$EQUATION, '<=>') #find the index (boolian)
equilibrum = filter(HMR_3_filter, equilibrum_index) #extract the equilbirum raw
non_equilbrium = filter(HMR_3_filter, !equilibrum_index)

#save data
save(equilibrum_index, equilibrum, non_equilbrium, file = "equilibrum_separated_dataframe.rda")


2. Among uni-directional group, separate left and right by '=>' using strsplit function
list = strsplit(non_equilbrium$EQUATION, '=>')
non_equ_dataFrame = data.frame(matrix = matrix(unlist(list), nrow=length(list), ncol=2, byrow=TRUE)) #dataframe

#comment : understandung the data type of input and output of function is important. (e.g.) vector, list, data.frame, matrix


3. Get boolian result table (left group, right group)
Boolian_ATP_consuming = str_detect(non_equ_dataFrame$consuming, "ATP") #list
Boolian_ATP_producing = str_detect(non_equ_dataFrame$producing, "ATP") #list

#Does there is only one ATP?
check = Boolian_ATP_consuming + Boolian_ATP_producing
unique(check) #NO
index_2 = which(check == 2) 

2#overlapping : 36 is ATP consuming. The others are all conversion
41 463 464 465 466 467 468 469 470 471 472 473 474 475 476 477 478 479 485 486 489
Boolian_ATP_producing[36] = FALSE
Boolian_ATP_consuming[index_2] = FALSE
Boolian_ATP_producing[index_2] = FALSE

4. get factor table. 1 stands for ATP-consuming. 2 stands fir ATP-producing

#make a new vector
name_factor = function(consuming_vector, producing_vector){
    #Assume : length of consuming_vector and producing_vector are same.
    #data type : all vector
    length = length(consuming_vector)
    print(length)
    result_vector = rep(NA, length)
    for (i in 1:length){
        #print(i)
        if (consuming_vector[i] == TRUE & producing_vector[i] == FALSE){
            result_vector[i] = 1 #consuming
         } else if (consuming_vector[i] == FALSE & producing_vector[i] == TRUE){
            result_vector[i] = 2 #producing
        } else if (consuming_vector[i] == FALSE & producing_vector[i] == FALSE){
            result_vector[i] = 0
         #conversion
        } else {
            result_vector[i] = 3 #check
        }
    }
    #print(i)
    return(result_vector)
}

factor_vector = name_factor(Boolian_ATP_consuming, Boolian_ATP_producing)
unique(factor_vector)

#Make a new version of data frame
non_equ_dataFrame_verion2 = mutate(non_equ_dataFrame, non_equilbrium$RXNID, non_equilbrium$`EC-NUMBER`, non_equilbrium$`GENE ASSOCIATION`,non_equilbrium$COMPARTMENT, non_equilbrium$SUBSYSTEM)
rm(non_equ_dataFrame)
non_equ_dataFrame_verion3 = mutate(non_equ_dataFrame_verion2, factor_vector)
non_equ_dataFrame_verion4 = non _equ_dataFrame_verion3[,c(3,1,2,4,5,6,8,7)] #change the sequence of data frame
colnames(non_equ_dataFrame_verion4) = c('RXNID', 'comsuming', 'producing', 'gene', 'compartment', 'subsystem', 'EC_number', 'factor')
#remove unnecessary memory space
rm(non_equ_dataFrame_verion2)
rm(factor_vector)
rm(Boolian_ATP_consuming)
rm(Boolian_ATP_producing)
rm(name_factor)

#save(equilibrum, non_equ_dataFrame_verion4, non_equilbrium, file = "equilibrum_and_non-equibrium_dataFrame.rda")
save(equilibrum, non_equ_dataFrame_verion4, file = "equilibrum_and_non-equibrium_dataFrame_2.rda")






<Extract gene list per each group>
#consuming group
geneList_consuming = non_equ_dataFrame_verion4 %>% 
    filter(factor == 1) %>%
    select(gene)
geneList_consuming = drop_na(geneList_consuming)
geneList_consuming = strsplit(geneList_consuming$`gene`, ';')  #convert dataframe into list 
geneList_consuming = unlist(geneList_consuming) # convert list into vector
save(geneList_consuming, file = "geneList_consuming.rda")
write(geneList_consuming, "geneList_consuming.txt")    

#producing group
geneList_producing = non_equ_dataFrame_verion4 %>% 
    filter(factor == 2) %>%
    select(gene)
geneList_producing = drop_na(geneList_producing)
geneList_producing = strsplit(geneList_producing$`gene`, ';')  #convert dataframe into list 
geneList_producing = unlist(geneList_producing) # convert list into vector
save(geneList_producing, file = "geneList_producing.rda")
write(geneList_producing, "geneList_producing.txt")





<sort RXN into equilbirum, ATP consuming, ATP producing>
#load data
load(file = "HMR_DataFrame.rda")
HMR_filter = HMR %>% drop_na('GENE ASSOCIATION')

#sort
1. equilibrum vs non-equilibrium
equilibrum_index = str_detect(HMR_filter$EQUATION, '<=>')
length(equilibrum_index[equilibrum_index == TRUE])
HMR_filter$COMMENTS[equilibrum_index] = 0

2. Among uni-directional group, separate left and right by '=>' using strsplit function
non_equilbrium_df = filter(HMR_filter, !equilibrum_index)
list = strsplit(non_equilbrium_df$EQUATION, '=>')
non_equ_dataFrame = data.frame(matrix = matrix(unlist(list), nrow=length(list), ncol=2, byrow=TRUE)) #dataframe
colnames(non_equ_dataFrame) = c('consuming', 'producing')
rm(list)

#comment : understandung the data type of input and output of function is important. (e.g.) vector, list, data.frame, matrix

3. Get boolian result table (left group, right group)
Boolian_ATP_consuming = str_detect(non_equ_dataFrame$consuming, "ATP") #list
Boolian_ATP_producing = str_detect(non_equ_dataFrame$producing, "ATP") #list

#Does there is only one ATP?
check = Boolian_ATP_consuming + Boolian_ATP_producing
unique(check) #NO
index_2 = which(check == 2) 
rm(check, index_2) # if no overlapping

#overlapping : 36 is ATP consuming. The others are all conversion
41 463 464 465 466 467 468 469 470 471 472 473 474 475 476 477 478 479 485 486 489
181  204 4256 4257 4258 4262 4263 4264 4265 4270 4273 4277
4280 4286 4293 4295 4300 4301 4310

Boolian_ATP_consuming[index_2] = FALSE
Boolian_ATP_producing[index_2] = FALSE
Boolian_ATP_consuming[181] = TRUE

#assign the sorted value into HMR dataframe
unique(non_equilbrium_df$COMMENTS)
non_equilbrium_df$COMMENTS[Boolian_ATP_consuming] = 1
non_equilbrium_df$COMMENTS[Boolian_ATP_producing] = 2
non_equilbrium_df = drop_na(non_equilbrium_df, COMMENTS)
rm(Boolian_ATP_consuming, Boolian_ATP_producing, equilibrum_index, non_equ_dataFrame)
save(non_equilbrium_df, file='non_equilibrum_df.rda')
save(HMR_filter, file='HMR_filter.rda')
    #length(non_equilbrium_df$COMPARTMENT[non_equilbrium_df$COMMENTS==1]) #374
    #length(non_equilbrium_df$COMPARTMENT[non_equilbrium_df$COMMENTS==2]) #5

<to make a map from gene to RXNID>
gene_to_RXNID_df = select(non_equilbrium_df, 'RXNID', 'GENE ASSOCIATION', 'COMMENTS')
colnames(gene_to_RXNID_df) = c('RXNID', 'gene', 'comment')
#convert_df_2 = separate(data = convert_df, col = gene, sep=';', into=LETTERS[1:30])
gene_list = strsplit(gene_to_RXNID_df$gene, ';')  #convert dataframe into list 

convert_gene_into_RXNID = function(list, HMR){
    #result_df = data.frame(gene=NA, RXNID=NA)
    gene_vector = c(NA)
    RXNID_vector = c(NA)
    factor = c(NA)
    for (i in 1:length(list)){
        for (j in 1:length(list[[i]])){
            #add_row(result_df,gene=i, RXNID=j)
            #add_row(result_df,gene=list[[i]][j], RXNID=conversion_df$'RXNID'[i])
            gene_vector = c(gene_vector, list[[i]][j])
            RXNID_vector = c(RXNID_vector, HMR$'RXNID'[i])
            factor= c(factor,HMR$'comment'[i])
            #print(i, j)
        }
    }
    result_df = data.frame(gene=gene_vector[2:length(gene_vector)], RXNID=RXNID_vector[2:length(gene_vector)], factor=factor[2:length(gene_vector)])
    return(result_df)
}
gene_to_RXNID_result = convert_gene_into_RXNID(gene_list, gene_to_RXNID_df)

save(gene_to_RXNID_df, gene_to_RXNID_result, non_equilbrium_df, file='gene_to_RXNID_result.rda')
rm(gene_list, gene_to_RXNID_df,non_equilbrium_df, convert_gene_into_RXNID)

#gene to HMR data
gene_to_HMRdata_df = left_join(x = gene_to_RXNID_result, y = HMR_3_filter, by = c('RXNID' = 'RXNID'))
save(gene_to_HMRdata_df, file='gene_to_HMRdata_df.rda')
