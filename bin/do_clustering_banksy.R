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
# add_pie_charts <- function(plot, proportions, cell_type_colors) {
#     # Scale proportions$size to c(4, 15)
#     proportions$size <- scales::rescale(proportions$size, to = c(0.3, 0.6))
#     plot <- plot + 
#         geom_scatterpie(aes(x, y, r=size),
#                         data = proportions, 
#                         cols = names(cell_type_colors),  # Use the dynamically created cell types
#                         pie_scale = 0.1) + 
#         scale_fill_manual(values = cell_type_colors) + # Use the dynamically created colors
#         coord_fixed()

#     return(plot)
# }

npcs <- 25
npcs_small <- 25

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
    norm_matrix <- na.omit(norm_matrix)
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

    seu <- RunBanksy(seu, lambda = 0.2, verbose = TRUE, assay = "RNA", slot = "data", features = 'variable', k_geom = 15, dimx = 'x', dimy = 'y')
    seu <- RunPCA(seu, assay="BANKSY", npcs = npcs, features = hvg)
    seu <- FindNeighbors(seu, dims = 1:npcs)

    for (res in c(0.1, 0.25, 0.5, 0.75, 1, 1.2)) {

        seu <- FindClusters(seu, resolution = res)
        seu[[paste0('banksy_', res)]] <- seu[['seurat_clusters']]

    }

    # proportions_banksy <- seu@meta.data %>%
    #     pivot_longer(cols = starts_with("banksy_"), names_to = "resolution", values_to = "cluster") %>%
    #     group_by(resolution, cluster, CellType) %>%
    #     summarise(Count = n()) %>%
    #     group_by(resolution, cluster) %>%
    #     mutate(Proportion = Count/sum(Count)) %>%
    #     ungroup()

    # # Get unique cell types and create a color palette
    # cell_types <- unique(proportions_banksy$CellType)
    # n_cell_types <- length(cell_types)

    # # Generate colors using the RColorBrewer palette (or any other palette)
    # colors <- brewer.pal(min(n_cell_types, 12), "Set3") # 'Set3' works well with categorical data
    # if (n_cell_types > 12) {
    #     colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_cell_types)
    # }

    # # Create a named vector for the colors
    # cell_type_colors <- setNames(colors, cell_types)

    # seu_aux <- seu
    # seu_aux[[]] <- na.omit(seu_aux[[]])
    # gbank <- clustree(seu_aux, 'banksy_')

    # # # Extract node positions (x, y) from the clustree plot
    # node_positions <- gbank$data[, c("x", "y", "banksy_", "cluster", "node", "size")] %>%
    #     rename(resolution = banksy_) %>%
    #     # Add louvain_ to the cluster name to avoid conflicts with other clustering methods
    #     mutate(resolution = paste0("banksy_", resolution))

    # # Reshape the data from long to wide format
    # proportions_wide_banksy <- proportions_banksy %>%
    #     pivot_wider(names_from = CellType,    # CellType becomes the new column names
    #                 values_from = Proportion, # The values will be from the 'Proportion' column
    #                 values_fill = list(Proportion = 0)) %>% # Fill missing proportions with 0
    #     group_by(resolution, cluster) %>%
    #     summarise(across(everything(), sum)) %>% # Sum the proportions for each cluster
    #     ungroup() %>%
    #     left_join(node_positions, by = c("resolution", "cluster")) # Merge with the node positions


    # gbank <- add_pie_charts(gbank, proportions_wide_banksy, cell_type_colors)
    # ggsave(paste0(output_plot_path, '/clustree_banksy_', top_, '_HVG_', method, '_', npcs, 'pcs' ,'.png'), height = 12, width = 20, limitsize = FALSE)

    # Save metadata with clustering results
    write.table(seu[[]], 
                paste0(output_csv_path, '/', 'top', top_, '_HVG_', method, '_', npcs,'pcs', '_clusters_banksy.csv'), 
                sep = ',', row.names = T, col.names = NA, quote = F)


}
