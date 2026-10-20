## Builds `data/HPA_ref.rda`: the single-cell type atlas of the Human Protein
## Atlas (average expression, in nCPM, of each gene in each cell type), as a
## SummarizedExperiment that SingleR can use as a reference. The values are
## stored in an assay named `logcounts`, the assay SingleR looks for. Run from
## the root of the repository.

library(dplyr)
library(tidyr)
library(tibble)
library(SummarizedExperiment)

zip <- tempfile(fileext = ".zip")
download.file("https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip", destfile = zip)
HPA_dat <- read.delim(unz(zip, "rna_single_cell_type.tsv"))[, -1] |>
    as_tibble() |>
    mutate(nCPM = as.numeric(nCPM)) |>
    pivot_wider(
        names_from  = Cell.type,
        values_from = nCPM,
        values_fn   = mean,
        values_fill = 0
    ) |>
    column_to_rownames("Gene.name")
HPA_cnts <- as.matrix(HPA_dat)
HPA_cnts[is.na(HPA_cnts)] <- 0

HPA_ref <- SummarizedExperiment(
    assays = list(logcounts = HPA_cnts),
    colData = DataFrame(label = colnames(HPA_dat))
)
save(HPA_ref, file = "data/HPA_ref.rda", compress = "xz")
