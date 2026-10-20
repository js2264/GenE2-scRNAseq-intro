# syntax=docker/dockerfile:1
ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
COPY . /opt/pkg

## Keep R's user cache (the book's `micromamba` and its conda env) inside the image
ENV R_USER_CACHE_DIR=/opt/R-cache

## Optionally pin `quarto`, e.g. to a release that renders `llms.txt` for books
## (quarto >= 1.11). Empty: keep the quarto shipped with the Bioconductor image.
## The download is retried on transient errors (e.g. HTTP 503 from GitHub)
ARG QUARTO_VERSION=""
ARG TARGETARCH
RUN if [ -n "${QUARTO_VERSION}" ]; then \
      curl -fsSL --retry 5 --retry-delay 10 -o /tmp/quarto.tar.gz "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH:-amd64}.tar.gz" && \
      mkdir -p "/opt/quarto-${QUARTO_VERSION}" && \
      tar -xzf /tmp/quarto.tar.gz -C "/opt/quarto-${QUARTO_VERSION}" --strip-components=1 && \
      ln -sf "/opt/quarto-${QUARTO_VERSION}/bin/quarto" /usr/local/bin/quarto && \
      rm /tmp/quarto.tar.gz ; \
    fi && \
    quarto --version

## Install book package. GITHUB_PAT arrives as a BuildKit secret, so that 
## GitHub-only dependencies do not hit the anonymous API rate limit, and is 
## never stored in the image
RUN --mount=type=secret,id=GITHUB_PAT \
    GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT 2>/dev/null || true)" \
    Rscript -e 'install.packages("remotes") ; repos <- BiocManager::repositories() ; remotes::install_local(path = "/opt/pkg/", repos=repos, dependencies=TRUE, build_vignettes=FALSE, upgrade=TRUE) ; sessioninfo::session_info(installed.packages()[,"Package"], include_base = TRUE)'

## Install the book's python environment, if any page runs python, and point
## every R session in the image at it
RUN if grep -rqsE '^```\{python' /opt/pkg/inst; then \
      Rscript -e 'BiocBook::setup_python(requirements = "/opt/pkg/inst/requirements.yml"); cat(sprintf("RETICULATE_PYTHON=%s\n", reticulate::py_config()$python), file = file.path(R.home("etc"), "Renviron.site"), append = TRUE)' && \
      chmod -R a+rX "${R_USER_CACHE_DIR}" && \
      Rscript -e 'stopifnot(nzchar(Sys.getenv("RETICULATE_PYTHON"))); cat("book python:", reticulate::py_config()$python, "\n")' ; \
    fi

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg

## ---- Workshop: make the image ready for RStudio Server (e.g. on Orchestra) ----

## Pre-download the data, so that nobody waits for a download on the day: the
## count matrix (already cached by the book build above) and the celldex
## reference used in the exercises. Both live in R's user cache, under
## R_USER_CACHE_DIR, which the image's init scripts pass on to RStudio
## sessions: make it writable by the RStudio user, and also reachable from
## ~/.cache/R, for R sessions started without R_USER_CACHE_DIR
RUN Rscript -e 'BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(ask = FALSE), "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112013/suppl/GSE112013_Combined_UMI_table.txt.gz") ; invisible(celldex::HumanPrimaryCellAtlasData())' && \
    chmod -R a+rwX "${R_USER_CACHE_DIR}" && \
    mkdir -p /home/rstudio/.cache && \
    ln -sfn "${R_USER_CACHE_DIR}/R" /home/rstudio/.cache/R && \
    chown -h rstudio:rstudio /home/rstudio/.cache /home/rstudio/.cache/R

## A copy of the workshop page and its figures in the RStudio user's home
RUN mkdir -p /home/rstudio/workshop && \
    cp /opt/pkg/inst/pages/workshop.qmd /home/rstudio/workshop/ && \
    cp -r /opt/pkg/inst/pages/figures /home/rstudio/workshop/ && \
    chown -R rstudio:rstudio /home/rstudio/workshop

## Orchestra keeps each participant's work in /data, and resets the rest of the
## container at every session: a start-up script makes /data writable by the
## RStudio user, and ~/data points to it, in RStudio's Files pane
RUN install -m 755 /opt/pkg/docker/cont-init.d/10_workshop_data /etc/cont-init.d/10_workshop_data && \
    ln -sfn /data /home/rstudio/data && \
    chown -h rstudio:rstudio /home/rstudio/data
