options(timeout = 600)

ensure_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

ensure_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

ensure_github <- function(repo, pkg = sub(".*/", "", repo)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(repo, dependencies = TRUE, repos = BiocManager::repositories())
  }
}

cran_packages <- c("remotes", "data.table", "readxl", "dplyr", "ggplot2")
invisible(lapply(cran_packages, ensure_cran))

ensure_cran("BiocManager")
BiocManager::install(version = "3.21", ask = FALSE)

bioc_packages <- c(
  "GenomicFeatures",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "GenomicRanges"
)
invisible(lapply(bioc_packages, ensure_bioc))

ensure_github("Townsend-Lab-Yale/cancereffectsizeR", pkg = "cancereffectsizeR")
ensure_github("Townsend-Lab-Yale/ces.refset.hg38@*release", pkg = "ces.refset.hg38")

message("Dependency installation finished.")
