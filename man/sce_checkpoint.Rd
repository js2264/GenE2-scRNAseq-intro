\name{sce_checkpoint}
\alias{sce_checkpoint}
\docType{data}
\title{
The workshop's SingleCellExperiment at the break
}
\description{
The \code{SingleCellExperiment} of the hands-on session as it is at the break,
for participants who fell behind: the adult human testis dataset (Guo et al.
2018, GEO accession GSE112013) after quality control, normalisation, feature
selection, PCA, UMAP and batch correction across donors with MNN.
}
\usage{data("sce_checkpoint")}
\format{
A \code{SingleCellExperiment} with the \code{counts} assay (to keep the
object small, the \code{logcounts} are not stored: compute them again with
\code{scrapper::normalizeRnaCounts.se()}, see the examples); the cell metadata (\code{Donor}, \code{Barcode}, \code{Replicate},
\code{Sample}), the quality metrics and the size factors in \code{colData};
the variance of each gene and the highly variable genes (\code{hvg}) in
\code{rowData}; and the \code{PCA}, \code{UMAP}, \code{MNN} and
\code{UMAP_MNN} reduced dimensions.
}
\source{
Guo et al. (2018), GEO accession
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112013}, processed
by \code{inst/scripts/sce_checkpoint.R}.
}
\references{
Guo J, et al. (2018). The adult human testis transcriptional cell atlas.
\emph{Cell Research} 28:1141-1157.
}
\examples{
data(sce_checkpoint)
sce_checkpoint
if (requireNamespace("scrapper", quietly = TRUE)) {
    sce <- scrapper::normalizeRnaCounts.se(
        sce_checkpoint, size.factors = sce_checkpoint$sum, block = sce_checkpoint$Sample
    )
    sce
}
}
\keyword{datasets}
