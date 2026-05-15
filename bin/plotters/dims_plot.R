library(Seurat)
library(ggplot2)


# Read terminal parameters
args <- commandArgs(trailingOnly = TRUE)
metadata_matrix <- args[1]
dims_matrix <- args[2]
output_plot_path <- args[3]
annotations_paths <- args[4:length(args)]

# Read metadata matrix
meta <- read.csv(metadata_matrix, row.names = 1)
dims <- read.csv(dims_matrix, row.names = 1)

# Filename 
filename <- strsplit(dims_matrix, split = "/")[[1]][length(strsplit(dims_matrix, split = "/")[[1]])]
top_ <- strsplit(filename, split = "_")[[1]][1]
norm_method <- strsplit(filename, split = "_")[[1]][3]
dim_method <- strsplit(filename, split = "_")[[1]][5]
dim_method <- substring(dim_method, 1, nchar(dim_method)-4)

# Read the correct annotations
for (annotations_path in annotations_paths) {
    if (grepl(paste0(top_, '_HVG_', norm_method), annotations_path)) {
        annotations <- read.csv(annotations_path, row.names = 1)
        break
    }
}

meta$CellType <- annotations$sup.clust

# Create a dummy matrix for seurat
dummy <- matrix(0, nrow = 5, ncol = nrow(dims))
colnames(dummy) <- rownames(dims)
seu <- CreateSeuratObject(counts = dummy, project = 'Seurat', assay="RNA", meta.data = meta)
seu[[dim_method]] <- CreateDimReducObject(embeddings = as.matrix(dims), assay = "RNA")

# Plot
p1 <- DimPlot(seu, reduction = dim_method) + NoLegend()
p2 <- FeaturePlot(seu, reduction = dim_method, features = 'Area')
p3 <- FeaturePlot(seu, reduction = dim_method, features = 'total_counts')
p4 <- FeaturePlot(seu, reduction = dim_method, features = 'fov')
# If the metadata has a column named donor_block_id, plot by donor_block_id
if ('donor_block_id' %in% colnames(meta)) {
    p5 <- DimPlot(seu, reduction = dim_method, group.by = 'donor_block_id')
} else {
    p5 <- DimPlot(seu, reduction = dim_method, group.by = 'CellType')
}
p6 <- DimPlot(seu, reduction = dim_method, group.by = 'CellType')

p <- CombinePlots(plots = list(p1, p2, p3, p4, p5, p6), ncol = 3)
ggsave(paste0(output_plot_path, '/', dim_method, '/',dim_method, '_top_', top_, '_HVG_', norm_method,'.png'), width = 30, height = 10)