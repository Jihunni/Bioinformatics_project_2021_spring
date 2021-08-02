library(tidyverse)
HCC = gsa.HCC_noHepatitis.summary
GTEx = gsa.GTEx_noHepatitis.summary

colnames(HCC)
HCC_merge = select(HCC, Name, 'p adj (dist.dir.up)', 'p (dist.dir.dn)')
colnames(HCC_merge) = c('Name', 'adjP_up_HCC', 'adjP_down_HCC')
GTEx_merge = select(GTEx, Name, 'p adj (dist.dir.up)', 'p (dist.dir.dn)')
colnames(GTEx_merge) = c('Name', 'adjP_up_GTEx', 'adjP_down_GTEx')

merge = inner_join(HCC_merge, GTEx_merge)
save(merge,file='sex_GTEx_merge_2021.05.15')
write.table(merge, file = "sex_GTEx_merge_2021.05.15.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE, append=FALSE)

rm(GTEx, HCC, merge, GTEx_merge, HCC_merge, input_reporterMetabolite_adjP, input_reporterMetabolite_lfc,sex)
