\name{HPA_ref}
\alias{HPA_ref}
\docType{data}
\title{
Human Protein Atlas single-cell type reference
}
\description{
A \code{SummarizedExperiment} holding the average expression (nCPM) of each
gene in 154 cell types, from the single-cell type atlas of the Human Protein
Atlas. The workshop uses it as a reference to annotate the clusters of the
adult human testis dataset with \code{SingleR}.
}
\usage{data("HPA_ref")}
\format{
A \code{SummarizedExperiment} with one assay, \code{logcounts} (genes x cell
types, the name of the assay \code{SingleR} looks for; the values are nCPM),
and one \code{colData} column, \code{label}, giving the name of each cell
type.
}
\source{
\url{https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip},
processed by \code{inst/scripts/HPA_ref.R}.
}
\references{
Karlsson M, et al. (2021). A single-cell type transcriptomics map of human
tissues. \emph{Science Advances} 7:eabh2169.
}
\examples{
data(HPA_ref)
HPA_ref
head(HPA_ref$label)
}
\keyword{datasets}
