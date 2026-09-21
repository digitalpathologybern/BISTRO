#!/usr/bin/env Rscript
# BISTRO -- pinned R environment installer
#
# envs/install_R_packages.R installs DIRECT dependencies and lets CRAN and Bioconductor
# resolve everything else against whatever is current, so it does not reproduce the
# environment the manuscript results were computed in. This file does.
#
# Source of truth: envs/R_packages_installed.csv, generated on 2026-09-10 from the library
# that produced the published results
# (/storage/research/igmp_dp_workspace/carreno_jose/R/x86_64-pc-linux-gnu-library/4.4),
# 371 packages, R R version 4.4.2 (2024-10-31), Bioconductor 3.20.
#
# renv is not available on this cluster, so the pin is driven from the manifest instead of
# from a renv.lock. Bioconductor packages are pinned by pinning the Bioconductor RELEASE,
# which is how Bioconductor versioning works; CRAN packages are pinned per package version.
#
# Usage:  Rscript envs/install_R_packages_pinned.R [--check]
#         --check reports what is missing or at the wrong version and installs nothing.

BIOC_VERSION <- "3.20"
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args

here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) "envs")
man_path <- file.path(here, "R_packages_installed.csv")
if (!file.exists(man_path)) stop("manifest not found: ", man_path)
man <- read.csv(man_path, stringsAsFactors = FALSE)

ip <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
ip$lib_rank <- match(ip$LibPath, .libPaths())
ip <- ip[order(ip$Package, ip$lib_rank), ]
ip <- ip[!duplicated(ip$Package), ]           # first .libPaths() hit is what R actually loads
have <- setNames(ip$Version, ip$Package)

status <- data.frame(
  package = man$package,
  want    = man$version,
  got     = unname(ifelse(man$package %in% names(have), have[man$package], NA)),
  stringsAsFactors = FALSE)
status$state <- ifelse(is.na(status$got), "MISSING",
                ifelse(status$got == status$want, "ok", "VERSION MISMATCH"))

cat(sprintf("Bioconductor release pinned to %s\n", BIOC_VERSION))
cat(sprintf("manifest: %d packages\n", nrow(status)))
print(table(status$state))
bad <- status[status$state != "ok", ]
if (nrow(bad)) print(head(bad[, c("package", "want", "got", "state")], 40), row.names = FALSE)

if (check_only) {
  quit(status = if (nrow(bad)) 1 else 0)
}
if (!nrow(bad)) { cat("environment already matches the manifest\n"); quit(status = 0) }

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")

bioc_pkgs <- rownames(available.packages(repos = BiocManager::repositories()[["BioCsoft"]]))
for (i in seq_len(nrow(bad))) {
  p <- bad$package[i]; v <- bad$want[i]
  if (p %in% bioc_pkgs) {
    message(sprintf("[bioc %s] %s -> %s", BIOC_VERSION, p, v))
    BiocManager::install(p, ask = FALSE, update = FALSE)
  } else {
    message(sprintf("[cran] %s -> %s", p, v))
    try(remotes::install_version(p, version = v, repos = "https://cloud.r-project.org", upgrade = "never"))
  }
}
cat("done; re-run with --check to verify\n")
