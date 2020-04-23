## -----------------------------------------------------------------------------
library("artdeco")
library("GEOquery")

# obtain path to testdata
gse73338 = getGEO('GSE73338')
transcriptome_data = exprs(gse73338$GSE73338_series_matrix.txt.gz)
gene_names = featureData(gse73338$GSE73338_series_matrix.txt.gz)[[3]]
rownames(transcriptome_data) = gene_names

# testrun
deconvolution_results = Deconvolve_transcriptome(
    transcriptome_data = transcriptome_data,
    deconvolution_algorithm = "nmf"
)

# show results
deconvolution_results[1:2,]

