#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scran)
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
output_csv_path   <- args[2]

set.seed(42)

# -----------------------------------------------
# 1. Stream-load counts efficiently
# -----------------------------------------------
message("Loading ", count_matrix_path)

dt <- fread(count_matrix_path, check.names = FALSE)
cell_names <- dt[[1]]
dt[[1]] <- NULL

for (j in seq_len(ncol(dt))) set(dt, j = j, value = as.numeric(dt[[j]]))

counts <- Matrix(as.matrix(dt), sparse = TRUE)
rownames(counts) <- as.character(cell_names)
rm(dt); gc()

counts <- Matrix::t(counts)
gc()
message("Counts loaded: ", nrow(counts), " genes x ", ncol(counts), " cells")

# -----------------------------------------------
# 2. Build SCE object (sparse)
# -----------------------------------------------
sce <- SingleCellExperiment(assays = list(counts = counts))
rm(counts); gc()
message("SCE built")

# -----------------------------------------------
# 3. Subsample for size factor estimation
# -----------------------------------------------
n_cells <- ncol(sce)
N_SUB   <- min(50000, n_cells)

if (n_cells > N_SUB) {
  message("Subsampling ", N_SUB, " / ", n_cells, " cells for size factor estimation")

  idx     <- sample(n_cells, N_SUB)
  sce_sub <- sce[, idx]

  clusters_sub <- quickCluster(
    sce_sub,
    block        = 5000,
    block.BPPARAM = MulticoreParam(workers = 5)
  )
  message("Clusters computed: ", length(unique(clusters_sub)))

  sce_sub <- computeSumFactors(sce_sub, clusters = clusters_sub, min.mean = 0.1)
  sf_sub  <- sizeFactors(sce_sub)
  message("Subset size factors computed (range: ",
          format(min(sf_sub)), " - ", format(max(sf_sub)), ")")

  # Extrapolate to full dataset via library-size ratio
  lib_sub <- colSums(assay(sce_sub, "counts"))
  lib_all <- colSums(assay(sce, "counts"))
  ratio   <- median(sf_sub / lib_sub)
  sf_all  <- lib_all * ratio

  sizeFactors(sce) <- sf_all
  rm(sce_sub, clusters_sub, sf_sub, lib_sub, idx); gc()

} else {
  message("Dataset small enough, running full deconvolution")

  clusters <- quickCluster(
    sce,
    block        = 5000,
    block.BPPARAM = MulticoreParam(workers = 5)
  )
  sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1)
  rm(clusters); gc()
}

sf <- sizeFactors(sce)
message("Final size factors (range: ",
        format(min(sf)), " - ", format(max(sf)), ")")

# -----------------------------------------------
# 4. Sparse normalization
# -----------------------------------------------
message("Normalizing counts with scran size factors (sparse)...")
norm_counts <- Matrix::t(Matrix::t(assay(sce, "counts")) / sf)

# -----------------------------------------------
# 5. Save size factors
# -----------------------------------------------
filename <- sub("\\.csv$", "", basename(count_matrix_path))

sf_path <- file.path(output_csv_path, paste0(filename, "_scranSF.csv"))
fwrite(data.table(cell = colnames(sce), size_factor = sf), sf_path)
message("Size factors saved: ", sf_path)

rm(sce); gc()

# -----------------------------------------------
# 6. Chunked CSV write (avoids full dense matrix)
# -----------------------------------------------
out_path   <- file.path(output_csv_path, paste0(filename, "_scran.csv"))
gene_names <- rownames(norm_counts)
block_size <- 5000

message("Writing normalized counts in chunks of ", block_size, " cells...")

# Write header row
writeLines(paste(c("cell", gene_names), collapse = ","), out_path)

n_total  <- ncol(norm_counts)
n_blocks <- ceiling(n_total / block_size)

for (b in seq_len(n_blocks)) {
  i_start <- (b - 1) * block_size + 1
  i_end   <- min(b * block_size, n_total)

  block <- as.matrix(norm_counts[, i_start:i_end, drop = FALSE])
  block_dt <- as.data.table(t(block), keep.rownames = TRUE)

  fwrite(block_dt, out_path, append = TRUE, col.names = FALSE)

  rm(block, block_dt); gc()

  if (b %% 10 == 0 || b == n_blocks) {
    message("  Written ", i_end, " / ", n_total, " cells")
  }
}

rm(norm_counts); gc()
message("Normalized counts saved: ", out_path)
message("Done.")