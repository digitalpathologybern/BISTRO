suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(future)
  library(Matrix)
  library(data.table)
})

options(future.globals.maxSize = 30000 * 1024^2) # Increase max size for future plan

# -------- Parse input arguments --------
args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
output_csv_path <- args[2]

message("Loading count matrix...")
dt <- fread(count_matrix_path)
cellnames <- dt[[1]]
dt[[1]] <- NULL
counts <- as.matrix(dt)
rownames(counts) <- as.character(cellnames)
counts <- Matrix(t(counts), sparse = TRUE)  # Transpose to genes x cells
rm(dt); gc()
message("Counts loaded (genes x cells): ", paste(dim(counts), collapse = " x "))
message("Min gene count in the dataset: ", min(colSums(counts)))


# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts))
message("SingleCellExperiment object created.")

seu <- as.Seurat(sce, counts = 'counts', data = 'counts')
message("Seurat object created.")

seu <- SCTransform(seu, vst.flavor = "v2", verbose = TRUE, assay = 'originalexp', variable.features.n = 20000)
message("SCTransform normalization completed.")

# Get input name and add _scTransform 
filename <- strsplit(count_matrix_path, split = "/")[[1]][length(strsplit(count_matrix_path, split = "/")[[1]])]
filename <- strsplit(toString(filename), split='.csv')

write.table(round(t(as.matrix(GetAssayData(object = seu, layer = "scale.data"))), 5), 
            paste0(output_csv_path, '/', filename, '_scTransform.csv'), 
            sep = ',', row.names = T, col.names = T, quote = F)