set.seed(123)   
library(scGSVA)   
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(parallel)
library(reshape2)
# module load CMake NEEDED


custom_msigdbr <- function(species_name, category_name, subcategory_name) {
    # Read csv
    pathways_path <- '/storage/homefs/jc23k573/~/cosmx-tma-data/code/cosmx-metastasis-tma/cosmx-analysis-files/gene-sets/celltypes_markers_v2024.1.csv'
    pathways <- read.csv(pathways_path)
    #   print(head(pathways))
    set_table <- pathways %>%
        rename('genes' = 'GeneID', 'category_code' = 'gs_cat', 'sub_category_code' = 'gs_subcat', 'set_name' = 'Annot') %>%
        # rename('gene_symbol' = 'genes', 'gs_cat' = 'category_code', 'gs_subcat' = 'sub_category_code', 'gs_name' = 'set_name') %>%
        # rename('GeneID' = 'genes', 'gs_cat' = 'category_code', 'gs_subcat' = 'sub_category_code', 'Annot' = 'set_name') %>%
        filter(organism == species_name & gs_cat == category_name & gs_subcat == subcategory_name)

    # Reorganize the table so we have one row per gene with all the other information
    set_table <- set_table %>%
        separate_rows(GeneID, sep = ',') %>%
        mutate(GeneID = trimws(GeneID)) %>%
        filter(GeneID != "")

    # Add one column with the number of genes per gs_cat
    set_table <- set_table %>%
        group_by(Annot) %>%
        mutate(total_n_genes = n()) %>%
        ungroup()

    return(set_table)
  
}


# Read terminal parameters
args <- commandArgs(trailingOnly = TRUE)
hvg_csv <- args[1]
output_plot_path <- args[2]
all_expression_csv <- args[3:length(args)]

# Read genes csv
hvg_df <- read.csv(hvg_csv, row.names = 1)
gene_id <- rownames(hvg_df)

# Get pathways
hsko_custom <- custom_msigdbr(species="Homo sapiens", category = "C4", subcategory = '3CA')
hsko_custom <- data.frame(hsko_custom) %>%
               select(GeneID, gs_cat, Annot)

top_ <- strsplit(hvg_csv, split = "/")[[1]][length(strsplit(hvg_csv, split = "/")[[1]])]
top_ <- strsplit(top_, split = "_")[[1]][2]

all_results <- data.frame(matrix(ncol = length(as.list(unique(hsko_custom['Annot']))$Annot), nrow = 0))
colnames(all_results) <- as.list(unique(hsko_custom['Annot']))$Annot

all_methods <- c()

for (expression_csv in all_expression_csv) {
    # expression_csv <- all_expression_csv[3]
    method <- strsplit(expression_csv, split = "/")[[1]][length(strsplit(expression_csv, split = "/")[[1]])]
    method <- strsplit(method, split = "_")[[1]][length(strsplit(method, split = "_")[[1]])]
    method <- substring(method, 1, nchar(method)-4)

    all_methods <- append(all_methods, method)

    # Read expression matrix
    data <- read.csv(expression_csv, row.names = 1)

    # Find the index name of True elements based on column method in hvg_df
    hvg_list <- rownames(hvg_df[as.logical(hvg_df[,method]), ])

    data <- data[,hvg_list]
    data <- t(as.matrix(data))
    
    # data <- data[,1:20000]
    # Pseudobulk per core
    core_id <- colnames(data)
    ## For each core_id, remove everything unitl _ and get the last part as an integer
    core_id <- sapply(core_id, function(x) as.numeric(strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])]))
    colnames(data) <- core_id
    ## Average expression per core
    pseudobulk_data <- aggregate(. ~ core_id, data.frame(t(data), core_id), mean)
    ## Transpose the result back so that genes are rows and coreIDs are columns
    pseudobulk_data <- t(pseudobulk_data[, -1])
    ## Add gene names back
    rownames(pseudobulk_data) <- rownames(data)
    data <- pseudobulk_data
    print(dim(data))

    # This function does not work properly for float data and sce or seu object, it expects count data or a data.frame
    res <- scgsva(data, hsko_custom, batch = 10000, ssgsea.norm = FALSE, parallel.sz = 0)
    res_df <- data.frame(res)
    
    # Calculate the mean value of the GSVA scores for each method
    method_mean <- colMeans(res_df)
    df_mean <- data.frame(t(method_mean))
    rownames(df_mean) <- method

    # Ensure subset_df has the same columns as all_results, with NA for missing columns
    subset_df_full <- data.frame(matrix(NA, nrow = 1, ncol = length(colnames(all_results))))
    colnames(subset_df_full) <- colnames(all_results)
    rownames(subset_df_full) <- method
    
    # Copy available data from subset_df to subset_df_full
    existing_cols <- intersect(colnames(df_mean), colnames(subset_df_full))
    subset_df_full[method, existing_cols] <- df_mean[1, existing_cols]
    
    # Combine with the all_results
    all_results <- rbind(all_results, subset_df_full)
}

all_results$Method <- rownames(all_results)
df_long <- melt(all_results, id.vars = "Method", variable.name = "Pathway", value.name = "MeanRepresentation")

n_methods <- length(all_methods)

g <- ggplot(df_long, aes(x = MeanRepresentation, y = Pathway, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  
  # Add the central vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +
  
  # Add horizontal lines for each pathway
  geom_hline(aes(yintercept = as.numeric(Pathway) - 0.17*n_methods), color = "gray80", linetype = "solid") +

  # Set plot labels and theme
  labs(title = "Pathway Representation by Normalization Method",
       x = "Mean Pathway Representation",
       y = "Pathway") +
  theme_minimal(base_size = 15) +
  
  # Customize colors and background
  scale_fill_brewer(palette = "Set1") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Remove the grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Increase the margin for the y-axis labels (pathway names)
    # axis.text.y = element_text(margin = margin(b = 5)),  # Adjust the margin for more space

    # Increase the plot margins (space around the entire plot)
    # plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave(paste0(output_plot_path, '/GSVA_analysis_top_', top_, 'HVG_CellTypes.png'), width = 15, height = 80, limitsize = FALSE)






