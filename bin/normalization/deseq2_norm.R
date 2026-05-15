#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(Matrix)
  library(data.table)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
output_csv_path <- args[2]

# ---- Load efficiently ----
message("Loading count matrix: ", count_matrix_path)
dt <- fread(count_matrix_path, nThread = 4, showProgress = TRUE)
cellnames <- dt[[1]]
dt[[1]] <- NULL

# Convert to sparse matrix (dgCMatrix)
counts <- as.matrix(dt)
rownames(counts) <- as.character(cellnames)
counts <- t(counts)
# counts <- Matrix(as.matrix(dt), sparse = TRUE)

rm(dt); gc()

message("Counts loaded: ", paste(dim(counts), collapse = " x "))

# ---- Estimate DESeq2 size factors ----
# DESeq2 requires dense matrix → use small subset if too large
if (object.size(counts) > 40e9) {
  message("⚠ Large matrix detected: converting temporarily to dense for DESeq2 normalization.")
}
size_factors <- DESeq2::estimateSizeFactorsForMatrix(counts, type = "poscounts")
names(size_factors) <- colnames(counts)

# ---- Normalize counts ----
message("Normalizing counts with DESeq2 size factors...")
norm_counts <- sweep(counts, 2, size_factors, FUN = "/")
rm(counts); gc()

# ---- Save outputs ----
filename <- sub("\\.csv$", "", basename(count_matrix_path))

# Save normalized counts in MatrixMarket format (fast & compact)
# out_norm <- file.path(output_csv_path, paste0(filename, "_deseq2.mtx"))
# writeMM(norm_counts, out_norm)
norm_dense <- as.matrix(norm_counts)
fwrite(as.data.table(t(norm_dense), keep.rownames = TRUE),
        file.path(output_csv_path, paste0(filename, "_deseq2.csv")))
message("Normalized counts saved.")

# Save size factors
sf_path <- file.path(output_csv_path, paste0(filename, "_deseq2SF.csv"))
fwrite(data.table(cell = names(size_factors), size_factor = size_factors), sf_path)

message("DESeq2 normalization complete ✅")

