library(InSituType)
library(SingleR)
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(Matrix)
library(data.table)

average_expression_by_celltype <- function(seurat_object, prefix) { 
  # Calculate average expression for each cell type
  seu.aggregated <- AggregateExpression(seurat_object, return.seurat = TRUE, group.by = "CellType")
  avg_expr_df <- GetAssayData(seu.aggregated, assay='RNA', layer='data')
  
  # Convert to data.frame and set the gene names as rownames
  # avg_expr_df <- as.data.frame(avg_expr)

  # Add prefix to the column names to indicate dataset source
  colnames(avg_expr_df) <- paste0(prefix, colnames(avg_expr_df))
  
  return(avg_expr_df)
}

clean_gene_names <- function(names_vec) {
  # Replace "-", "_", " ", ":", "/" with "."
  return(gsub("[-_:/ ]", ".", names_vec))
}

# Read terminal parameters
args <- commandArgs(trailingOnly = TRUE)
reference_dataset_path <- args[1]
count_matrix_path <- args[2]
metadata_matrix_path <- args[3]
hvg_csv_path <- args[4]
output_csv_path <- args[5]
norm_matrices_path <- args[6:length(args)]

# Load reference data
# Check if the reference dataset is an RDS file or a RData object
if (grepl(".rds", reference_dataset_path)) {
  ref.expr <- readRDS(reference_dataset_path)
  # Prepare the dataset used for annotation
  ref.avg <- average_expression_by_celltype(ref.expr, "Reference_")
} else {
  load(reference_dataset_path)
  ref.avg <- as.matrix(profile_matrix)
}

rownames(ref.avg) <- clean_gene_names(rownames(ref.avg))

# Read data to be annotated
# count.expr <- t(read.csv(count_matrix_path, row.names = 1))
dt <- fread(count_matrix_path)
cellnames <- dt[[1]]
dt[[1]] <- NULL
count.expr <- Matrix(as.matrix(dt), sparse = TRUE)
rownames(count.expr) <- as.character(cellnames)
# count.expr <- t(count.expr)  # Transpose to genes x cells
rm(dt); gc()

meta <- fread(metadata_matrix_path)
cellnames_meta <- meta[[1]]
meta[[1]] <- NULL
norm.meta <- as.data.frame(meta)
rownames(norm.meta) <- as.character(cellnames_meta)
rm(meta); gc()

colnames(count.expr) <- clean_gene_names(colnames(count.expr))
# norm.meta <- read.csv(metadata_matrix_path, row.names = 1)

# Read hvg csv
hvg_df <- read.csv(hvg_csv_path, row.names = 1)
rownames(hvg_df) <- clean_gene_names(rownames(hvg_df))

top_ <- strsplit(hvg_csv_path, split = "/")[[1]][length(strsplit(hvg_csv_path, split = "/")[[1]])]
top_ <- strsplit(top_, split = "_")[[1]][2]

# unique_cell_types <- c()

for (norm_matrix_path in norm_matrices_path) {
  filename <- strsplit(norm_matrix_path, split = "/")[[1]][length(strsplit(norm_matrix_path, split = "/")[[1]])]
  filename <- strsplit(filename, split = "_")[[1]][length(strsplit(filename, split = "_")[[1]])]
  method <- substr(filename, start = 1, stop = nchar(filename)-4)

  print(method)
  print(filename)

  method <- gsub("-", ".", method)

  # Prepare CosMx for the annotation
  # Count matrix
  # count.expr.t <- t(count.expr)

  hvg_list <- rownames(hvg_df[as.logical(hvg_df[,method]), ])
  
  # Verify intersection to prevent crashes
  common_genes <- intersect(hvg_list, colnames(count.expr))
  common_genes <- intersect(common_genes, rownames(ref.avg))
  
  if (length(common_genes) == 0) {
    stop("No common genes found between Reference, Counts, and HVG list after cleaning names!")
  }

  # If special characters in gene names of count.expr, replace them with "."
  # if (any(grepl("\\.", colnames(count.expr)))) {
  # if (any(grepl("-", colnames(count.expr)))) {
  #   colnames(count.expr) <- gsub("-", ".", colnames(count.expr))
  # } else if (any(grepl('_', colnames(count.expr)))) {
  #   colnames(count.expr) <- gsub("_", ".", colnames(count.expr))
  # } else if (any(grepl(" ", colnames(count.expr)))) {
  #   colnames(count.expr) <- gsub(" ", ".", colnames(count.expr))
  # } else if (any(grepl(":", colnames(count.expr)))) {
  #   colnames(count.expr) <- gsub(":", ".", colnames(count.expr))
  # } else if (any(grepl("/", colnames(count.expr)))) {
  #   colnames(count.expr) <- gsub("/", ".", colnames(count.expr))
  # }

  # if (any(grepl("-", hvg_list))) {
  #   hvg_list <- gsub("-", ".", hvg_list)
  # } else if (any(grepl('_', hvg_list))) {
  #   hvg_list <- gsub("_", ".", hvg_list)
  # } else if (any(grepl(" ", hvg_list))) {
  #   hvg_list <- gsub(" ", ".", hvg_list)
  # } else if (any(grepl(":", hvg_list))) {
  #   hvg_list <- gsub(":", ".", hvg_list)
  # } else if (any(grepl("/", hvg_list))) {
  #   hvg_list <- gsub("/", ".", hvg_list)
  #  }
  
  print(hvg_list)
  print(colnames(count.expr))
  print(common_genes)
  
  print(paste("Using", length(common_genes), "genes for annotation."))

  count.expr.small <- count.expr[,common_genes]
  print('Matrix subsetted successfully')

  # Vector of mean neg probes per cell
  negmean <- norm.meta$total_counts_Negative/20
  names(negmean) <- rownames(norm.meta)

  sup <- insitutypeML(
    x = count.expr.small,
    neg = negmean, 
    reference_profiles = ref.avg[common_genes,]
  )

  # Annotate the CosMx dataset
  annotation <- data.frame(sup$clust)

  write.table(annotation, 
              paste0(output_csv_path, '/', top_, '_HVG_', method, '_annotation.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)

  # unique_cell_types <- unique(c(unique_cell_types, unique(sup$clust)))
}

# unique_cell_types <- unique(unique_cell_types)
# color_list <- read.csv('/storage/homefs/jc23k573/~/cosmx-tma-data/code/cosmx-metastasis-tma/color_list.csv', sep = "\t", stringsAsFactors = FALSE, header = FALSE)
# integers <- seq(from=1, to=length(color_list[,1]), by=1)
# random_colors <- sample(integers, length(unique_cell_types))
# # random_colors <- sample(color_list[,1], length(unique_cell_types))

# color_palette <- data.frame(
#   cell_type = unique_cell_types,
#   color = random_colors
# )

# # Save the color palette to a CSV file
# write.csv(color_palette, 'color_palette.csv', row.names = FALSE)
