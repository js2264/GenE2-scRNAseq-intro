<!-- badges: start -->
📦 [Repo](https://github.com/js2264/GenE2-scRNAseq-intro) [![rcmdcheck](https://img.shields.io/github/actions/workflow/status/js2264/GenE2-scRNAseq-intro/rcmdcheck.yml?label=Package%20check)](https://github.com/js2264/GenE2-scRNAseq-intro/actions/workflows/rcmdcheck.yml)   
📖 [Book](https://js2264.github.io/GenE2-scRNAseq-intro/devel) [![deployment](https://img.shields.io/github/actions/workflow/status/js2264/GenE2-scRNAseq-intro/pages/pages-build-deployment?label=Book%20deployment)](https://github.com/js2264/GenE2-scRNAseq-intro/actions/workflows/pages/pages-build-deployment)  
🐳 [Docker](https://github.com/js2264/GenE2-scRNAseq-intro/pkgs/container/gene2scrnaseqintro) [![biocbook](https://img.shields.io/github/actions/workflow/status/js2264/GenE2-scRNAseq-intro/biocbook.yml?label=Docker%20image)](https://github.com/js2264/GenE2-scRNAseq-intro/actions/workflows/biocbook.yml)  
<!-- badges: end -->

# Introduction to single-cell RNA-seq analysis with Bioconductor

**A 2h30 live-coding workshop for the GenE2 Master 2, *Single-Cell Omics and
Spatial Transcriptomics* teaching unit.**

📖 **Book: <https://js2264.github.io/GenE2-scRNAseq-intro/devel/>**

This repository is a [`BiocBook`](https://github.com/js2264/BiocBook): an R
package (`GenE2scRNAseqIntro`) whose `inst/` folder holds a Quarto book. Every
push to `devel` builds a Docker image with everything the workshop needs,
renders the book inside it, and publishes both.

| Chapter | Content | Format |
|---|---|---|
| Welcome | Goals, and how to reach the workshop environment | read |
| 1. From cells to counts | How a 10x Genomics experiment turns a tissue into a UMI count matrix | read before the session, ~15 min |
| 2. Hands-on | Count matrix → `SingleCellExperiment` → QC → normalisation → HVGs → PCA/UMAP → batch correction → clustering → markers → `SingleR` annotation → pseudobulk `DESeq2`, on one page | live, 2h30 |

The dataset is the adult human testis atlas of
[Guo *et al.* (2018)](https://doi.org/10.1038/s41422-018-0099-2), GEO
[GSE112013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112013):
6,490 cells, 3 donors × 2 technical replicates.

## Running the material

**On Orchestra**, during the session: one RStudio Server per participant, from
the image below.

**With Docker**: everything is pre-installed, and the data are already
downloaded.

```sh
docker run -e PASSWORD=gene2 -p 8787:8787 ghcr.io/js2264/gene2scrnaseqintro:devel
```

then open <http://localhost:8787> and log in as `rstudio` / `gene2`. A copy of
the workshop page is in `~/workshop/`.

**Locally**, with R >= 4.6.0, a recent Bioconductor and about 8 GB of RAM:

```r
BiocManager::install("js2264/GenE2-scRNAseq-intro", dependencies = TRUE)
```

## For the instructor

### How the book is built

| Push to | GitHub Actions does |
|---|---|
| `devel` | `biocbook`: builds the Docker image against Bioconductor devel, renders the book inside it, pushes `ghcr.io/js2264/gene2scrnaseqintro:devel` (also tagged `latest` and with the package version) and deploys the book to `gh-pages`, under `docs/devel/`. `rcmdcheck`: `R CMD check` of the package. |
| a pull request to `devel` | `rcmdcheck` only. |
| any other branch, by hand (*Actions* > *biocbook* > *Run workflow*) | `biocbook` builds the image and renders the book, without pushing or deploying anything. Use it to check a branch before merging. |

Edit the book in `inst/`: the pages are `inst/index.qmd` and `inst/pages/*.qmd`,
their order is set in `inst/assets/_book.yml`. To preview it locally:

```r
bb <- BiocBook::BiocBook(".")
BiocBook::preview(bb)
```

### One-off settings, before the first deployment

1. **Pages**: *Settings* > *Pages* > *Build and deployment*: deploy from the
   branch `gh-pages`, folder `/docs`. The book is then served at
   <https://js2264.github.io/GenE2-scRNAseq-intro/devel/>.
2. **Make the Docker image public.** GHCR may create the package as private,
   and `docker pull` (Orchestra's included) then fails without credentials.
   After the first successful build, go to
   <https://github.com/users/js2264/packages/container/gene2scrnaseqintro/settings>
   > *Danger Zone* > *Change visibility* > **Public**, and check with a
   `docker pull` from a logged-out machine.

### Orchestra

Point the workshop at the image `ghcr.io/js2264/gene2scrnaseqintro:devel` (or
at the tag of a given package version, to freeze the environment for the day),
with RStudio on port `8787`. Set `DISABLE_AUTH=true` as an environment variable
to skip the RStudio login. Reading the dense count table needs about 4 GB of
RAM at its peak, so pick a size with at least 8 GB.

The image is ready for participants:

- the count matrix and the `celldex` reference of the exercises are already
  cached, so nothing is downloaded during the session;
- `~/workshop/` holds a copy of the workshop page and its figures.

### Data shipped with the package

| Object | What | Made by |
|---|---|---|
| `HPA_ref` | Human Protein Atlas single-cell-type reference, for `SingleR` | `inst/scripts/HPA_ref.R` |
| `sce_checkpoint` | The `SingleCellExperiment` at the break, for participants who fall behind | `inst/scripts/sce_checkpoint.R` |

Re-run the scripts from the root of the repository if the workflow changes.

## Licence and citation

MIT (see `LICENSE`). If you reuse this workshop, please keep the attribution
and cite the source of the data:

> Guo J, Grow EJ, Mlcochova H, Maher GJ, Lindskog C, *et al.* (2018).
> The adult human testis transcriptional cell atlas.
> *Cell Research* **28**:1141–1157.

Adapted from the
[BiocAfrica 2026 scRNA-seq workshop](https://github.com/js2264/BiocAfrica_scRNAseq2026).
