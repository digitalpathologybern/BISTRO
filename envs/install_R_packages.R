#!/usr/bin/env Rscript
# ============================================================================
# BISTRO -- R dependency installer
# ----------------------------------------------------------------------------
# Installs the R packages the pipeline needs. Only DIRECT dependencies are
# listed here; CRAN and Bioconductor resolve the transitive ones. The full set
# of 371 packages present in the environment that produced the manuscript
# results is recorded separately in R_packages.csv, which is a provenance
# record rather than an install list.
#
# Usage:
#     Rscript envs/install_R_packages.R            # install what is missing
#     Rscript envs/install_R_packages.R --check    # report only, install nothing
#     Rscript envs/install_R_packages.R --optional # also install optional extras
#
# Requires R >= 4.4 and a working compiler toolchain; several Bioconductor
# packages build from source.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
with_optional <- "--optional" %in% args

# ---- Dependency sets -------------------------------------------------------
# Derived from the library() calls in the scripts that run_BISTRO.nf actually
# invokes. Version numbers are the ones used for the manuscript.

cran_pkgs <- c(
    "data.table",   # 1.18.2.1
    "dplyr",        # 1.2.1
    "future",       # 1.70.0
    "ggplot2",      # 4.0.3
    "Matrix",       # 1.7-5
    "msigdbr",      # 26.1.0
    "Seurat",       # 5.5.0
    "sctransform",  # 0.4.3
    "tibble",       # 3.3.1
    "tidyr",        # 1.3.2
    "XML",          # 3.99-0.23
    "xml2",         # 1.5.2
    "RColorBrewer", # 1.1-3
    "reshape2"      # 1.4.5
)

bioc_pkgs <- c(
    "BiocParallel",         # 1.40.2
    "DESeq2",               # 1.46.0
    "edgeR",                # 4.4.2
    "limma",                # 3.62.2
    "GOfuncR",              # 1.26.0
    "scran",                # 1.34.0
    "scuttle",              # 1.16.0
    "SingleCellExperiment", # 1.28.1
    "SingleR",              # 2.8.0
    "SpatialExperiment",    # 1.16.0
    "SpaNorm"               # 1.0.0
)

# Not on CRAN or Bioconductor. InSituType is required for the annotation step
# (skip_annotation = false, the default).
github_pkgs <- c(InSituType = "Nanostring-Biostats/InSituType")

# Only needed for the GSVA step, which is commented out in run_BISTRO.nf.
github_optional <- c(scGSVA = "guokai8/scGSVA")

# Auxiliary scripts in bin/ that the pipeline does NOT invoke
# (reduce_dimensions.R, do_clustering_banksy.R, do_clustering_lou_lei.R)
# additionally need Banksy, clustree, scatterpie, SeuratWrappers and
# SpatialPCA. They are deliberately not installed here.


# ---- Helpers ---------------------------------------------------------------

# "Installed" and "loadable" are different questions, and conflating them sends
# people down the wrong path. Several Bioconductor packages here link against
# system libraries -- igraph needs libglpk, magick needs libMagick++, sf and
# terra need GDAL and UDUNITS -- and when one of those is absent the R package
# is present on disk but its namespace will not load. Reporting that as
# "missing" would prompt a reinstall that cannot fix it.

is_installed <- function(p) length(find.package(p, quiet = TRUE)) > 0
is_loadable  <- function(p) requireNamespace(p, quietly = TRUE)

load_error <- function(p) {
    msg <- tryCatch({ loadNamespace(p); NULL },
                    error = function(e) conditionMessage(e))
    if (is.null(msg)) return(NA_character_)
    gsub("\\s+", " ", substr(msg, 1, 160))
}

# Returns only the genuinely absent packages; unloadable ones are reported
# separately because installing them again would not help.
report <- function(pkgs, label) {
    cat(sprintf("\n%s\n", label))
    missing <- character(0)
    broken <- character(0)
    for (p in pkgs) {
        if (!is_installed(p)) {
            cat(sprintf("  [ ] %-22s not installed\n", p))
            missing <- c(missing, p)
        } else if (!is_loadable(p)) {
            cat(sprintf("  [!] %-22s %s  INSTALLED BUT WILL NOT LOAD\n",
                        p, as.character(utils::packageVersion(p))))
            cat(sprintf("      %s\n", load_error(p)))
            broken <- c(broken, p)
        } else {
            cat(sprintf("  [x] %-22s %s\n", p,
                        as.character(utils::packageVersion(p))))
        }
    }
    list(missing = missing, broken = broken)
}


# ---- Check mode ------------------------------------------------------------

if (check_only) {
    cat("BISTRO R dependency check\n")
    cat(sprintf("R version: %s\n", getRversion()))

    r1 <- report(cran_pkgs, "CRAN")
    r2 <- report(bioc_pkgs, "Bioconductor")
    r3 <- report(names(github_pkgs), "GitHub (required)")
    invisible(report(names(github_optional), "GitHub (optional, GSVA only)"))

    missing_required <- c(r1$missing, r2$missing, r3$missing)
    broken_required <- c(r1$broken, r2$broken, r3$broken)

    cat("\n", strrep("=", 60), "\n", sep = "")
    if (length(missing_required) == 0 && length(broken_required) == 0) {
        cat("All required packages are installed and loadable.\n")
        quit(status = 0)
    }
    if (length(missing_required)) {
        cat(sprintf("NOT INSTALLED (%d): %s\n", length(missing_required),
                    paste(missing_required, collapse = ", ")))
        cat("  Run this script without --check to install them.\n")
    }
    if (length(broken_required)) {
        cat(sprintf("INSTALLED BUT NOT LOADABLE (%d): %s\n",
                    length(broken_required),
                    paste(broken_required, collapse = ", ")))
        cat("  These are a SYSTEM LIBRARY problem, not an R one -- reinstalling\n")
        cat("  the R package will not fix it. See the load errors above and the\n")
        cat("  system library table in the README's System requirements section.\n")
        cat("  On a module-based HPC the usual cause is simply not having loaded\n")
        cat("  the matching modules (GLPK, ImageMagick, GDAL, UDUNITS, libxml2).\n")
    }
    quit(status = 1)
}


# ---- Install ---------------------------------------------------------------

cat("BISTRO R dependency installation\n")
cat(sprintf("R version: %s\n", getRversion()))
cat(sprintf("Library:   %s\n\n", .libPaths()[1]))

if (getRversion() < "4.4") {
    warning("BISTRO was developed and tested on R 4.4.2. ",
            "Older versions are untested and Bioconductor 3.20 requires R 4.4.")
}

if (!is_installed("BiocManager")) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!is_installed("remotes")) {
    cat("Installing remotes...\n")
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

# BiocManager::install() handles CRAN and Bioconductor packages alike and keeps
# them consistent with the installed Bioconductor release.
to_install <- c(cran_pkgs, bioc_pkgs)
missing <- to_install[!vapply(to_install, is_installed, logical(1))]

if (length(missing)) {
    cat(sprintf("\nInstalling %d CRAN/Bioconductor package(s):\n  %s\n\n",
                length(missing), paste(missing, collapse = ", ")))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
} else {
    cat("\nAll CRAN/Bioconductor packages already present.\n")
}

install_gh <- function(pkgs) {
    for (nm in names(pkgs)) {
        if (is_installed(nm)) {
            cat(sprintf("  %s already installed\n", nm))
            next
        }
        cat(sprintf("  installing %s from %s ...\n", nm, pkgs[[nm]]))
        tryCatch(
            remotes::install_github(pkgs[[nm]], upgrade = "never"),
            error = function(e)
                message(sprintf("  FAILED to install %s: %s", nm, conditionMessage(e)))
        )
    }
}

cat("\nGitHub packages (required):\n")
install_gh(github_pkgs)

if (with_optional) {
    cat("\nGitHub packages (optional):\n")
    install_gh(github_optional)
} else {
    cat("\nSkipping optional packages (scGSVA). Pass --optional to include them.\n")
}

# ---- Summary ---------------------------------------------------------------

required <- c(cran_pkgs, bioc_pkgs, names(github_pkgs))
still_missing <- required[!vapply(required, is_installed, logical(1))]

cat("\n", strrep("=", 60), "\n", sep = "")
if (length(still_missing) == 0) {
    cat("All required R packages are installed.\n")
    cat("Verify at any time with:\n")
    cat("    Rscript envs/install_R_packages.R --check\n")
    quit(status = 0)
}

cat(sprintf("%d package(s) still missing: %s\n",
            length(still_missing), paste(still_missing, collapse = ", ")))
cat("Check the build log above for compilation errors; several Bioconductor\n")
cat("packages need system libraries (libxml2, GDAL, UDUNITS, GLPK).\n")
quit(status = 1)
