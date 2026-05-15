library(Banksy) # Banksy must be loaded before Seurat
library(Seurat)
library(SpatialPCA)
library(clustree)
library(SeuratWrappers)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scatterpie)

# Python is needed for leidenalg library
# To make this work, I had to modify the file: /storage/homefs/jc23k573/.local/easybuild/software/R/4.4.1-gfbf-2023b/lib64/R/etc/Renviron
# module load libxml2
# module load GLPK/5.0-GCCcore-10.3.0
# module load ImageMagick

# Save pca, umap, tsne coordinates + clustering trees + clustering results
# Function to create pie charts for each cluster
add_pie_charts <- function(plot, proportions, cell_type_colors) {
    # Scale proportions$size to c(4, 15)
    proportions$size <- scales::rescale(proportions$size, to = c(0.3, 0.6))
    plot <- plot + 
        geom_scatterpie(aes(x, y, r=size),
                        data = proportions, 
                        cols = names(cell_type_colors),  # Use the dynamically created cell types
                        pie_scale = 0.1) + 
        scale_fill_manual(values = cell_type_colors) + # Use the dynamically created colors
        coord_fixed()

    return(plot)
}

npcs <- 25
npcs_small <- 25
run_tsne <- FALSE

# Read terminal parameters
args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
metadata_matrix_path <- args[2]
hvg_csv_path <- args[3]
output_csv_path <- args[4]
output_plot_path <- args[5]
use_global <- args[6]
# norm_matrices_path <- args[6:length(args)]

# All remaining arguments are the normalized matrices or annotations
remaining_paths <- args[7:length(args)]

# Separate norm matrices and annotation files by recognizing naming patterns
norm_matrices_path <- remaining_paths[grepl("_filtered_filtered_counts_", remaining_paths)]  # Norm matrices have 'counts' in the name
annotation_files_path <- remaining_paths[grepl("_annotation", remaining_paths)]  # Annotation files have 'annotation' in the name

for (norm_matrix_path in norm_matrices_path) {

    filename <- strsplit(norm_matrix_path, split = "/")[[1]][length(strsplit(norm_matrix_path, split = "/")[[1]])]
    filename <- strsplit(filename, split = "_")[[1]][length(strsplit(filename, split = "_")[[1]])]
    method <- substr(filename, start = 1, stop = nchar(filename)-4)

    top_ <- strsplit(hvg_csv_path, split = "/")[[1]][length(strsplit(hvg_csv_path, split = "/")[[1]])]
    top_ <- strsplit(top_, split = "_")[[1]][2]

    # Load matrices
    count_matrix <- read.csv(count_matrix_path, header = TRUE, row.names = 1)
    norm_matrix <- read.csv(norm_matrix_path, header = TRUE, row.names = 1)
    hvg_csv <- read.csv(hvg_csv_path, header = TRUE, row.names = 1)
    meta <- read.csv(metadata_matrix_path, header = TRUE, row.names = 1)
    # Read the correct annotations
    for (annotations_path in annotation_files_path) {
        if (grepl(paste0(top_, '_HVG_', method), annotations_path)) {
            annotations <- read.csv(annotations_path, row.names = 1)
            break
        }
    }

    hvg <- rownames(hvg_csv[as.logical(hvg_csv[,method]),])

    # Create Seurat object
    seu <- CreateSeuratObject(counts = t(count_matrix), project = "cosmx", assay = "RNA")
    LayerData(seu, 'data') <- t(norm_matrix)
    VariableFeatures(seu) <- hvg
    if (use_global == '1'){
        seu[['x']] <- meta['CenterX_global_px']  + meta['fov_centerX'] / (0.12e-3)
        seu[['y']] <- meta['CenterY_global_px'] + meta['fov_centerY'] / (0.12e-3)
    } else {
        seu[['x']] <- meta['CenterX_global_px']
        seu[['y']] <- meta['CenterY_global_px']
    }
    seu[['CellType']] <- annotations['sup.clust']

    seu <- ScaleData(seu, features = VariableFeatures(object = seu))

    seu <- RunPCA(seu, 
                npcs = 50,
                features = VariableFeatures(object = seu))
    seu <- RunUMAP(seu, n.neighbors = 30, min.dist = 0.1, dims = 1:npcs, slot = "scale.data")

    if (run_tsne) {
        seu <- RunTSNE(seu, dims = 1:npcs_small, slot = "scale.data", check_duplicates = FALSE)
    }
    

    elbow <- ElbowPlot(seu, ndims = length(hvg))
    ggsave(paste0(output_plot_path, '/elbow_top_', top_, '_HVG_', method, '_', npcs, 'pcs' ,'.png'))

    # Get coordinates
    pca_coords <- Embeddings(seu, reduction = "pca")
    umap_coords <- Embeddings(seu, reduction = "umap")

    if (run_tsne) {
        tsne_coords <- Embeddings(seu, reduction = "tsne")
    }

    # Write to file
    write.table(pca_coords, 
                paste0(output_csv_path, '/pca/', top_, '_HVG_', method, '_', npcs,'pcs', '_pca.csv'), 
                sep = ',', row.names = T, col.names = NA, quote = F)
    write.table(umap_coords, 
                paste0(output_csv_path, '/umap/', top_, '_HVG_', method, '_', npcs,'pcs', '_umap.csv'), 
                sep = ',', row.names = T, col.names = NA, quote = F)
    if (run_tsne) {
        write.table(tsne_coords, 
                    paste0(output_csv_path, '/tsne/', top_, '_HVG_', method, '_', npcs,'pcs', '_tsne.csv'), 
                    sep = ',', row.names = T, col.names = NA, quote = F)
    }

}
