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

**On [Orchestra](https://orchestraplatform.org)**, during the session: one
RStudio Server per participant, from the image below (see
[Orchestra](#orchestra)).

**With Docker**: everything is pre-installed, and the data are already
downloaded.

```sh
docker run -v "$PWD":/data -e PASSWORD=gene2 -p 8787:8787 ghcr.io/js2264/gene2scrnaseqintro:devel
```

then open <http://localhost:8787> and log in as `rstudio` / `gene2`. A copy of
the workshop page is in `~/workshop/`, and the current folder of your computer
appears as `~/data`.

**Locally**, with R >= 4.6.0, a recent Bioconductor and about 8 GB of RAM:

```r
BiocManager::install("js2264/GenE2-scRNAseq-intro", dependencies = TRUE)
```

## For the instructor

### How the book is built

| Push to | GitHub Actions does |
|---|---|
| `devel` | `biocbook`: builds the Docker image against Bioconductor devel, renders the book inside it, pushes `ghcr.io/js2264/gene2scrnaseqintro:devel` (also tagged `latest`, with the package version, and with the commit's short SHA) and deploys the book to `gh-pages`, under `docs/devel/`. `rcmdcheck`: `R CMD check` of the package. |
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

The image meets Orchestra's
[no-auth contract](https://orchestraplatform.org/docs/user-guide/configuring-images/)
as it is: RStudio listens on port `8787`, and Orchestra sets `DISABLE_AUTH=true`
itself, so no environment variable or container argument is needed. For
participants, the image also:

- caches the count matrix and the `celldex` reference of the exercises, so
  that nothing is downloaded during the session;
- holds a copy of the workshop page and its figures in `~/workshop/`;
- makes `/data`, Orchestra's per-participant workspace, writable by the
  RStudio user at start-up, and links it as `~/data`. With a persistent
  workspace, what participants save there comes back in their next session;
  everything else is reset.

**Add the workshop to the catalog** with a
[*Workshop template* issue](https://github.com/orchestraplatform/orchestraplatform/issues/new?template=workshop-template.yml)
on the Orchestra repository:

| Field | Value |
|---|---|
| Display name | `Introduction to scRNA-seq analysis (GenE2 Master 2)` |
| Slug | `gene2-scrnaseq-intro` |
| Description | `2h30 live-coding introduction to single-cell RNA-seq analysis with Bioconductor, for the GenE2 Master 2 "Single-Cell Omics and Spatial Transcriptomics" teaching unit.` |
| Image | `ghcr.io/js2264/gene2scrnaseqintro:<short SHA>`: the 7-character SHA of the `devel` commit whose image you want |
| App port | `8787` |
| Size | Large (2 CPU, 8Gi memory) |
| Tags | bioconductor, rstudio |
| Environment variables, Container args | leave empty |
| Storage size | `10Gi` |
| Persistent workspace | Per-user (keep /data between sessions) |
| Landing URL | <https://js2264.github.io/GenE2-scRNAseq-intro/devel/> |
| Source repo URL | <https://github.com/js2264/GenE2-scRNAseq-intro> |

A bot validates the issue and comments; once a maintainer adds the `approve`
label, it opens a pull request with the template, and merging that pull
request deploys the workshop to the catalog.

- **Size.** A session peaks at about 4 GB of RAM, when it reads the dense
  count table: Standard (4Gi) is too tight, Large (8Gi) leaves room for the
  exercises.
- **Image tag.** Orchestra's nodes cache the images they pull, so a tag that
  is pushed again (`devel`, `latest`, or the package version until it is
  bumped) can give different participants different builds. The short-SHA tag
  of a commit never changes. To update the workshop, push to `devel`, then
  submit the issue form again with the same slug and the new short SHA. Sessions
  already running keep their image: update before the event, not during it.
- **On the day.** The default session length (4h) fits a 2h30 workshop. Tell
  the platform admin your headcount and start time ahead of the event, so that
  enough capacity is warmed up: the first launch on a cold node pulls the whole
  image, which takes a few minutes.

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
