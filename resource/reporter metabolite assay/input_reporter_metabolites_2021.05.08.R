BiocManager::install("piano", dependencies=TRUE)
library(piano)
?piano
?runGSA

GEM_liver_hepatocytes = loadGSC('liver_hepatocytes.xml')

mmrGsc´ë½Å hmrGsc

#HCC(male tumor w/o Hepatitis vs female tumor w/o Hepatitis)
gsa.HCC_noHepatitis = runGSA(input_reporterMetabolite_padj, input_reporterMetabolite_lfc, gsc=hmrGsc, geneSetStat="reporter", signifMethod="nullDist", nPerm=1000, gsSizeLim=c(5,100))
gsa.HCC_noHepatitis.summary = GSAsummaryTable(gsa.HCC_noHepatitis)
##network.HCC_noHepatitis <- networkPlot(gsa.HCC_noHepatitis, class="distinct", direction="both", significance=0.005, label="numbers")

save(gsa.HCC_noHepatitis, file='gsa.HCCnoH_GTEx_2021.05.15.rda')
save(gsa.HCC_noHepatitis.summary, file= 'gsa.gsa.HCCnoH_GTEx_2021.05.15.rda')
write.table(gsa.HCC_noHepatitis.summary, file = "gsa.HCCnoH_GTEx_summary_2021.05.08.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote=FALSE, append=FALSE)
