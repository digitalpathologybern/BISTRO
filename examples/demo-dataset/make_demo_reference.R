#!/usr/bin/env Rscript
# ============================================================================
# BISTRO demo dataset: build the InSituType reference
# ----------------------------------------------------------------------------
# Converts the genes x cell-types CSV written by make_demo_dataset.py into the
# .RData object that bin/annotation/run_insitutype.R expects.
#
# That script branches on the reference file extension:
#   *.rds    -> readRDS(), expects a Seurat object, averages by "CellType"
#   anything else -> load(), expects an object literally named `profile_matrix`
#
# The demo takes the second path, so the saved object MUST be named
# profile_matrix, and must be a matrix (or coercible by as.matrix) with gene
# symbols as row names and cell type names as column names.
#
# Usage:
#     Rscript make_demo_reference.R demo_reference_profiles.csv demo_reference.RData
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: Rscript make_demo_reference.R <profiles.csv> <output.RData>")
}

profiles_csv <- args[1]
output_rdata <- args[2]

message("Reading reference profiles: ", profiles_csv)

profile_matrix <- as.matrix(
    read.csv(profiles_csv, row.names = 1, check.names = FALSE)
)

message("  profile_matrix: ", nrow(profile_matrix), " genes x ",
        ncol(profile_matrix), " cell types")
message("  cell types: ", paste(colnames(profile_matrix), collapse = ", "))

save(profile_matrix, file = output_rdata)
message("Wrote ", output_rdata)
