## Builds `data/sce_checkpoint.rda`: the SingleCellExperiment of the hands-on
## session as it is at the break, for participants who fall behind. This is the
## code of `inst/pages/workshop.qmd` up to the break: keep both in sync. Run
## from the root of the repository.

library(BiocFileCache)
library(Matrix)
library(SingleCellExperiment)
library(scrapper)

bfc <- BiocFileCache(ask = FALSE)
umi_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112013/suppl/GSE112013_Combined_UMI_table.txt.gz"
umi <- read.delim(bfcrpath(bfc, umi_url), header = TRUE, row.names = 1)
sce <- SingleCellExperiment(assays = list(counts = Matrix(as.matrix(umi), sparse = TRUE)))
rm(umi)

split_names <- strsplit(colnames(sce), "\\.")
sce$Donor <- sapply(split_names, function(x) x[1])
sce$Barcode <- sapply(split_names, function(x) x[2])
sce$Replicate <- sapply(split_names, function(x) x[3])
sce$Sample <- paste(sce$Donor, sce$Replicate, sep = "_rep")

is_mito <- grepl("^MT-", rownames(sce))
sce <- quickRnaQc.se(sce, subsets = list(Mito = is_mito), block = sce$Sample)
sce <- sce[, sce$keep]
sce <- normalizeRnaCounts.se(sce, size.factors = sce$sum, block = sce$Sample)
sce <- chooseRnaHvgs.se(sce, block = sce$Sample, top = 2000)
sce <- runPca.se(sce, features = rowData(sce)$hvg, number = 25)
sce <- runUmap.se(sce)
sce <- correctMnn.se(sce, block = sce$Donor)
sce <- runUmap.se(sce, reddim.type = "MNN", output.name = "UMAP_MNN")

## The log-normalised values are computed on the fly from the counts and the
## size factors: storing them would store the counts a second time, so they
## are dropped, and re-computed after loading by the same call as above
sce_checkpoint <- sce
assay(sce_checkpoint, "logcounts") <- NULL
save(sce_checkpoint, file = "data/sce_checkpoint.rda", compress = "xz")
