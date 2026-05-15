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
    seu <- FindNeighbors(seu, dims = 1:npcs)

    # Perform clustering at different resolutions
    for (res in c(0.1, 0.25, 0.5, 0.75, 1, 1.2)) {

        seu <- FindClusters(seu, resolution = res, algorithm = 1)
        seu[[paste0('louvain_', res)]] <- seu[['seurat_clusters']]

        seu <- FindClusters(seu, resolution = res, algorithm = 4, method = 'igraph')
        seu[[paste0('leiden_', res)]] <- seu[['seurat_clusters']]
    }

    proportions_louvain <- seu@meta.data %>%
        pivot_longer(cols = starts_with("louvain_"), names_to = "resolution", values_to = "cluster") %>%
        group_by(resolution, cluster, CellType) %>%
        summarise(Count = n()) %>%
        group_by(resolution, cluster) %>%
        mutate(Proportion = Count/sum(Count))%>%
        ungroup()

    proportions_leiden <- seu@meta.data %>%
        pivot_longer(cols = starts_with("leiden"), names_to = "resolution", values_to = "cluster") %>%
        group_by(resolution, cluster, CellType) %>%
        summarise(Count = n()) %>%
        group_by(resolution, cluster) %>%
        mutate(Proportion = Count/sum(Count))%>%
        ungroup()

    # Get unique cell types and create a color palette
    cell_types <- unique(proportions_louvain$CellType)
    n_cell_types <- length(cell_types)

    # Generate colors using the RColorBrewer palette (or any other palette)
    colors <- brewer.pal(min(n_cell_types, 12), "Set3") # 'Set3' works well with categorical data
    if (n_cell_types > 12) {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_cell_types)
    }

    # Create a named vector for the colors
    cell_type_colors <- setNames(colors, cell_types)

    glou <- clustree(seu, 'louvain_', node_size_range = c(6, 20))

    # # Extract node positions (x, y) from the clustree plot
    node_positions <- glou$data[, c("x", "y", "louvain_", "cluster", "node", "size")] %>%
        rename(resolution = louvain_) %>%
        # Add louvain_ to the cluster name to avoid conflicts with other clustering methods
        mutate(resolution = paste0("louvain_", resolution))

    # Reshape the data from long to wide format
    proportions_wide_louvain <- proportions_louvain %>%
        pivot_wider(names_from = CellType,    # CellType becomes the new column names
                    values_from = Proportion, # The values will be from the 'Proportion' column
                    values_fill = list(Proportion = 0)) %>% # Fill missing proportions with 0
        group_by(resolution, cluster) %>%
        summarise(across(everything(), sum)) %>% # Sum the proportions for each cluster
        ungroup() %>%
        left_join(node_positions, by = c("resolution", "cluster")) # Merge with the node positions

    glou <- add_pie_charts(glou, proportions_wide_louvain, cell_type_colors)
    ggsave(paste0(output_plot_path, '/clustree_louvain_', top_, '_HVG_', method, '_', npcs, 'pcs' ,'.png'), height = 12, width = 20, limitsize = FALSE)

    ##########################################################
    # Get unique cell types and create a color palette
    cell_types <- unique(proportions_leiden$CellType)
    n_cell_types <- length(cell_types)

    # Generate colors using the RColorBrewer palette (or any other palette)
    colors <- brewer.pal(min(n_cell_types, 12), "Set3") # 'Set3' works well with categorical data
    if (n_cell_types > 12) {
    colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_cell_types)
    }

    # Create a named vector for the colors
    cell_type_colors <- setNames(colors, cell_types)

    glei <- clustree(seu, 'leiden_')

    # # Extract node positions (x, y) from the clustree plot
    node_positions <- glei$data[, c("x", "y", "leiden_", "cluster", "node", "size")] %>%
        rename(resolution = leiden_) %>%
        # Add louvain_ to the cluster name to avoid conflicts with other clustering methods
        mutate(resolution = paste0("leiden_", resolution))

    # Reshape the data from long to wide format
    proportions_wide_leiden <- proportions_leiden %>%
        pivot_wider(names_from = CellType,    # CellType becomes the new column names
                    values_from = Proportion, # The values will be from the 'Proportion' column
                    values_fill = list(Proportion = 0)) %>% # Fill missing proportions with 0
        group_by(resolution, cluster) %>%
        summarise(across(everything(), sum)) %>% # Sum the proportions for each cluster
        ungroup() %>%
        left_join(node_positions, by = c("resolution", "cluster")) # Merge with the node positions

    glei <- add_pie_charts(glei, proportions_wide_leiden, cell_type_colors)
    ggsave(paste0(output_plot_path, '/clustree_leiden_', top_, '_HVG_', method, '_', npcs, 'pcs' ,'.png'), height = 12, width = 20, limitsize = FALSE)

    # Save metadata with clustering results
    write.table(seu[[]], 
                paste0(output_csv_path, '/', 'top', top_, '_HVG_', method, '_', npcs,'pcs', '_clusters_loulei.csv'), 
                sep = ',', row.names = T, col.names = NA, quote = F)


}
