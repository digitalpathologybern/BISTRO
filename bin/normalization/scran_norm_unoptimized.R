#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scran)
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(magrittr)
  library(peakRAM)
})

args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
output_csv_path   <- args[2]

message("Loading ", count_matrix_path)

## ───────────────────────────────────────────────
## 1️⃣  Stream-load counts efficiently
## ───────────────────────────────────────────────
# fread is fast, but we keep everything as numeric right away
dt <- fread(count_matrix_path, check.names = FALSE)
cell_names <- dt[[1]]
dt[[1]] <- NULL

# Convert to numeric in-place (no copy)
for (j in seq_len(ncol(dt))) set(dt, j = j, value = as.numeric(dt[[j]]))

# Create a sparse matrix directly without transposing yet
counts <- Matrix(as.matrix(dt), sparse = TRUE)
rownames(counts) <- as.character(cell_names)
rm(dt); gc()

# Transpose lazily (avoid double memory)
counts <- Matrix::t(counts)
gc()
message("Counts loaded: ", nrow(counts), " genes × ", ncol(counts), " cells")

## ───────────────────────────────────────────────
## 3️⃣  Build SCE object (sparse)
## ───────────────────────────────────────────────
sce <- SingleCellExperiment(assays = list(counts = counts))
rm(counts); gc()
message("SCE loaded")

## ───────────────────────────────────────────────
## 4️⃣  Clustering with memory-aware parameters
## ───────────────────────────────────────────────
# quickCluster can use a block size (cells per block) to cap memory
# Larger blocks = faster but more memory
res <- peakRAM({
  clusters <- quickCluster(sce, block = 5000, block.BPPARAM = BiocParallel::MulticoreParam(workers = 5))
})
message("Peak RAM used during clustering: ", res)
message("Clusters calculated: ", length(unique(clusters)))
gc()

## ───────────────────────────────────────────────
## 5️⃣  Compute sum factors in chunks
## ───────────────────────────────────────────────
sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1)
message("Sum factors calculated")

sf <- sizeFactors(sce)
message("Size factors computed (range: ",
        format(min(sf)), "–", format(max(sf)), ")")


## ───────────────────────────────────────────────
## 6️⃣  Write output incrementally
## ───────────────────────────────────────────────
filename <- sub("\\.csv$", "", basename(count_matrix_path))

# Normalize in chunks to avoid full dense matrix in RAM
# out_path <- file.path(output_csv_path, paste0(filename, "_scran.csv"))
# con <- file(out_path, open = "wt")
# writeLines("cell," %>% paste(colnames(sce), collapse = ","), con)
# block_size <- 1000
# for (i in seq(1, 2, by = block_size)) {
# # for (i in seq(1, nrow(assay(sce)), by = block_size)) {
#   block <- assay(sce, "counts")[i:min(i + block_size - 1, nrow(sce)), , drop = FALSE]
#   norm_block <- t(t(as.matrix(block)) / sf)
#   write.table(norm_block, con, sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE, append = TRUE)
#   rm(block, norm_block); gc()
# }
# close(con)
message("Normalizing counts with scran size factors...")
norm_counts <- sweep(assay(sce, "counts"), 2, sf, FUN = "/")
rm(sce); gc()
# Save normalized counts
norm_dense <- as.matrix(norm_counts)
fwrite(as.data.table(t(norm_dense), keep.rownames = TRUE),
        file.path(output_csv_path, paste0(filename, "_scran.csv")))
message("Normalized counts saved.")

# Save size factors
sf_path <- file.path(output_csv_path, paste0(filename, "_scranSF.csv"))
fwrite(data.table(cell = names(sf), size_factor = sf), sf_path)
# data.table(cell = names(sf), size_factor = sf) |>
#   fwrite(file.path(output_csv_path, paste0(filename, "_scranSF.csv")))

message("Files saved successfully")