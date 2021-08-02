BiocManager::install("piano", dependencies=TRUE)
library(piano)
?piano
?runGSA

mmrGsc´ë½Å hmrGsc

#HCC(male tumor w/o Hepatitis vs female tumor w/o Hepatitis)
gsa.GTEx_noHepatitis = runGSA(input_reporterMetabolite_adjP, input_reporterMetabolite_lfc, gsc=hmrGsc, geneSetStat="reporter", signifMethod="nullDist", nPerm=1000, gsSizeLim=c(5,100))
gsa.GTEx_noHepatitis.summary = GSAsummaryTable(gsa.GTEx_noHepatitis)
##network.HCC_noHepatitis <- networkPlot(gsa.HCC_noHepatitis, class="distinct", direction="both", significance=0.005, label="numbers")

save(gsa.GTEx_noHepatitis, file='gsa.GTEx_2021.05.13.rda')
save(gsa.GTEx_noHepatitis.summary, file= 'gsa.HCC_GTEx_summary_2021.05.13.rda')
write.table(gsa.GTEx_noHepatitis.summary, file = "gsa.HCC_GTEx_summary_2021.05.13.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote=FALSE, append=FALSE)