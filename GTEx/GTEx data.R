#library
#install.packages('Cepa')
#library('CePa') #not updated
library(readr)
install.packages('phantasus')
BiocManager::install("phantasus")
library(phantasus)
#load GTEx data (gct format)
#GTEx = read.gct('./GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct')
#save(GTEx, file='GTEx_transcript_expected_count.rda')
#load('GTEx_gene_reads.rda')

# read_tsv(file = 'liverGTExCountTab.txt',
#     col_names = TRUE,
#     col_types = NULL,
#     locale = default_locale(),
#     na = c("", "NA"),
#     quoted_na = TRUE,
#     quote = "\"",
#     comment = "",
#     trim_ws = TRUE,
#     skip = 0,
#     n_max = Inf,
#     guess_max = min(1000, n_max),
#     progress = show_progress(),
#     skip_empty_rows = TRUE)

GTEx_liver = read_tsv("liverGTExCountTab.txt", header =TRUE, sep="")
save(GTEx_liver, file = 'GTEx_liver.rda')
#total 226 samples

dds <- DESeqDataSetFromMatrix(countData = gtex,
                              colData = sampledata,
                              design = ~ tissue) #generate the deseq data set

dds <- dds[ rowSums(counts(dds)) > 1, ] #remove genes with zero counts

vsd <- vst(dds, blind = FALSE) #normalization considering tissue
