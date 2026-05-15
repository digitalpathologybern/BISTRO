#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
output_csv_path <- args[2]
total_count <- as.numeric(args[3])

# Define letter suffix
letter <- switch(as.character(total_count),
                 "1e+06" = "m",
                 "10000" = "10k",
                 "100" = "100",
                 "custom")

message("Loading counts from: ", count_matrix_path)

# 1️⃣ Efficient loading
message("Loading count matrix...")
dt <- fread(count_matrix_path)
cellnames <- dt[[1]]
dt[[1]] <- NULL
counts <- Matrix(as.matrix(dt), sparse = TRUE)
rownames(counts) <- as.character(cellnames)
counts <- t(counts)  # Transpose to genes x cells
rm(dt); gc()
message("Counts loaded (genes x cells): ", paste(dim(counts), collapse = " x "))

# 2️⃣ Compute library sizes (total counts per cell)
# Columns = cells
lib_sizes <- colSums(counts)
sf_path <- file.path(output_csv_path, 
                     paste0(sub("\\.csv$", "", basename(count_matrix_path)), 
                            "_cp", letter, "SF.csv"))
fwrite(data.table(cell = names(lib_sizes), size_factor = lib_sizes), sf_path)
print(paste("Size factors saved to:", sf_path))

# 3️⃣ Compute CPM efficiently
# CPM = (counts / colSums(counts)) * total_count
# This works directly on sparse matrices:
cpm <- sweep(counts, 2, lib_sizes, FUN = "/") * as.numeric(total_count)
rm(counts); gc()

# 4️⃣ Save results in sparse MatrixMarket format (fast & compact)
filename <- sub("\\.csv$", "", basename(count_matrix_path))
# out_path <- file.path(output_csv_path, paste0(filename, "_cp", letter, ".mtx"))
# Matrix::writeMM(cpm, out_path)
cpm_dense <- as.matrix(cpm)
fwrite(as.data.table(t(cpm_dense), keep.rownames = TRUE),
        file.path(output_csv_path, paste0(filename, "_cp", letter, ".csv")))

message("Normalization completed successfully.")
