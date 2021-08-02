library('readr')
library(tidyverse)

load(GTEx_liver)

#GTEx dataframe (subjectID and sex) from liver transcriptome
GTEx_sample_id = data.frame(sample_id=colnames(GTEx_liver), conversion=colnames(GTEx_liver))
GTEx_sample_id$conversion = gsub("[.]","-", GTEx_sample_id$conversion)
GTEx_sample_id = cbind(GTEx_sample_id, subject_ID=substr
dim(GTEx_sample_id) #226 3

# #merge
# GTEx_sample_id = left_join(x=GTEx_sample_id, y=subject_sex, by = c('subject_ID' ='subjectId'))
# dim(GTEx_sample_id) #226 4
# table(GTEx_sample_id$sex) #male 24, female 44



###again with new data
liver_histology=read_tsv('GTEx_liver_histology.txt')
GTEx_sample_id = GTEx_sample_id[c('sample_id', 'conversion')]
GTEx_sample_id = cbind(GTEx_sample_id, liver_histology)
save(GTEx_sample_id, file='GTEx_liver_subject_sex.rda')
