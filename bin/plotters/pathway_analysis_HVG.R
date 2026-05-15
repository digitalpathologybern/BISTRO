library(GOfuncR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(msigdbr)
library(xml2)
library(XML)

# Read terminal parameters
args <- commandArgs(trailingOnly = TRUE)
genes_csv <- args[1]
output_plot_path <- args[2]
output_csv_path <- args[3]

# Read genes csv
filename <- strsplit(genes_csv, split = "/")[[1]][length(strsplit(genes_csv, split = "/")[[1]])]
filename <- strsplit(filename, split = "_")[[1]][2]
genes <- read.csv(genes_csv, row.names = 1)

# Bring the colnames to a column
genes <- data.frame(genes)
genes$gene_ids <- rownames(genes)

# Iterate through each column in genes except gene_ids and create input_hyper dataframe
# all_stats <- list()
# for (i in 1:(ncol(genes)-1)) {
#     method <- colnames(genes)[i]
#     input_hyper <- data.frame(genes$gene_ids, as.integer(as.logical(genes[,i])))
#     res_hyper <- go_enrich(input_hyper, n_randset = 1000)
#     stats <- res_hyper[[1]]
#     stats$method <- method

#     all_stats[[i]] <- stats
# }

# combined_stats <- bind_rows(all_stats)

# top_5_nodes_per_ontology <- combined_stats %>%
#   group_by(ontology, method) %>%
#   arrange(raw_p_overrep) %>%   # Sort by raw_p_overrep (smallest values first)
#   slice_head(n = 5) %>%        # Keep the top 5 nodes for each ontology
#   ungroup()

# top_5_nodes_per_ontology <- top_5_nodes_per_ontology %>%
#   mutate(node_name = factor(node_name, levels = unique(node_name[order(ontology)])))  # Add ontology to node_name

# # Create a plot using ggplot where the x axis is a method (represented in the columns of the genes dataframe) and the y axis are the different pathways, grouped by the ontology. Then plot a dot with the size being the inversed on raw_p_overrep, present in stats
# ggplot(top_5_nodes_per_ontology, aes(x = method, y = node_name, color = ontology, size = 1 - raw_p_overrep)) +
#   geom_point(alpha = 0.7) +  # Add points with some transparency
#   scale_size_continuous(name = "Inverse Raw p Overrep") +  # Legend title for size
#   scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +  # Add some space to the x-axis
#   theme_minimal() +  # A clean theme
#   labs(
#     x = "Normalization Method",
#     y = "Node Name",
#     title = "Dot Plot of Node Names by Normalization Method and Ontology"
#   ) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if necessary
#     plot.title = element_text(hjust = 0.5),              # Center the plot title
#     panel.background = element_rect(fill = "white", color = NA),  # White panel background
#     plot.background = element_rect(fill = "white", color = NA),   # White plot background
#     legend.background = element_rect(fill = "white", color = NA)  # White legend background
#   )

# # Save the plot and add filename
# ggsave(paste0(output_plot_path, '/GO_analysis.png'), width = 15, limitsize = FALSE)


#custom_msigdbr <- function(species_name, category_name, subcategory_name) {
#  # Read csv
#  pathways_path <- '/storage/homefs/jc23k573/main/cosmx-tma-data/code/cosmx-metastasis-tma/cosmx-analysis-files/gene-sets/celltypes_markers_v2024.1.csv'
#  pathways <- read.csv(pathways_path)
#
#  set_table <- pathways %>%
#    rename('gene_symbol' = 'genes', 'gs_cat' = 'category_code', 'gs_subcat' = 'sub_category_code', 'gs_name' = 'set_name') %>%
#    filter(organism == species_name & gs_cat == category_name & gs_subcat == subcategory_name)
#
#  # Reorganize the table so we have one row per gene with all the other information
#  set_table <- set_table %>%
#    separate_rows(gene_symbol, sep = ',') %>%
#    mutate(gene_symbol = trimws(gene_symbol)) %>%
#    filter(gene_symbol != "")
#
#  # Add one column with the number of genes per gs_cat
#  set_table <- set_table %>%
#    group_by(gs_name) %>%
#    mutate(total_n_genes = n()) %>%
#    ungroup()
#
#  # print(head(set_table))
#  return(set_table)
# 
#}

#######################
hallmarks <- msigdbr(species="Homo sapiens", collection = "H") #C7 - Sub: IMMUNESIGDB

hallmarks <- hallmarks %>%
    group_by(gs_name) %>%
    mutate(total_n_genes = n()) %>%
    ungroup()

results <- list()
for (i in 1:(ncol(genes)-1)){
    method <- colnames(genes)[i]
    input_genes <- genes$gene_ids[as.logical(genes[,i])]
    
    # Create a dataframe for this method
    method_data <- hallmarks %>%
      filter(gene_symbol %in% input_genes) %>%  # Filter for input genes
      distinct(gs_name, gene_symbol, .keep_all = TRUE) %>%
      group_by(gs_name, total_n_genes) %>%
      summarise(
        total_genes_overlap = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        percentage_overlap = (total_genes_overlap / total_n_genes) * 100,
        method = method
      ) %>%
      arrange(desc(percentage_overlap))
  
    # Add this method's data to the results list
    results[[method]] <- method_data
}

combined_stats <- bind_rows(results)

write.table(combined_stats, paste0(output_csv_path, '/pathway_analysis_top_', filename, '_HVG_Hallmarks.csv'), sep = ",", row.names = FALSE)

# Plot the data with ontologies on the y-axis and methods on the x-axis
ggplot(combined_stats, aes(x = method, y = gs_name, color = method, size = percentage_overlap)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Percentage of Overlap", range = c(2, 10)) + 
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +  # Adjust spacing if needed
  coord_cartesian(clip = "off") +  # Allow plot elements to extend beyond plot area
  labs(
    x = "Normalization Method",
    y = "Hallmark",
    title = "Percentage of detected HVGs in Hallmarks"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),             # Center the plot title
    plot.margin = margin(t = 20, r = 20, b = 50, l = 20),  # Adjust margins for space
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA),   # White plot background
    legend.background = element_rect(fill = "white", color = NA)  # White legend background
  )

ggsave(paste0(output_plot_path, '/Pathway_analysis_top_', filename, 'HVG_Hallmarks.png'), width = 15, height = 25, limitsize = FALSE)


#######################
# hallmarks <- custom_msigdbr(species="Homo sapiens", category = "C4", subcategory = '3CA') #C7 - Sub: IMMUNESIGDB

hallmarks <- msigdbr(species = "Homo sapiens", collection = "C4", subcollection = "3CA")

hallmarks <- hallmarks %>%
    group_by(gs_name) %>%
    mutate(total_n_genes = n()) %>%
    ungroup()

results <- list()
for (i in 1:(ncol(genes)-1)){
    method <- colnames(genes)[i]
    input_genes <- genes$gene_ids[as.logical(genes[,i])]
    
    # Create a dataframe for this method
    method_data <- hallmarks %>%
      filter(gene_symbol %in% input_genes) %>%
      distinct(gs_name, gene_symbol, .keep_all = TRUE) %>%
      group_by(gs_name, total_n_genes) %>%
      summarise(
        total_genes_overlap = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        percentage_overlap = (total_genes_overlap / total_n_genes) * 100,
        method = method
      ) %>%
      arrange(desc(percentage_overlap))
  
    # Add this method's data to the results list
    results[[method]] <- method_data
}


combined_stats <- bind_rows(results)

write.table(combined_stats, paste0(output_csv_path, '/pathway_analysis_top_', filename, '_HVG_CellAtlas.csv'), sep = ",", row.names = FALSE)

# Plot the data with ontologies on the y-axis and methods on the x-axis
ggplot(combined_stats, aes(x = method, y = gs_name, color = method, size = percentage_overlap)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Percentage of Overlap", range = c(2, 10)) + 
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +  # Adjust spacing if needed
  coord_cartesian(clip = "off") +  # Allow plot elements to extend beyond plot area
  labs(
    x = "Normalization Method",
    y = "Cell type pathway",
    title = "Percentage of detected HVGs in Cell types pathways"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5),             # Center the plot title
    plot.margin = margin(t = 20, r = 20, b = 50, l = 20),  # Adjust margins for space
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA),   # White plot background
    legend.background = element_rect(fill = "white", color = NA)  # White legend background
  )
# Save the plot and add filename

ggsave(paste0(output_plot_path, '/Pathway_analysis_top_', filename, 'HVG_CellAtlas.png'), width = 15, height = 25, limitsize = FALSE)

############################ I DONT GET INTERESTING RESULTS WITH THIS ############################
######## Plot the percentage of overlap of my HVG wiht the different immune pathways and create a 3 columns plot

# hallmarks <- msigdbr(species="Homo sapiens", category = "C7",subcategory="IMMUNESIGDB") #H
# results <- list()
# for (i in 1:(ncol(genes)-1)){
#     method <- colnames(genes)[i]
#     input_genes <- genes$gene_ids[as.logical(genes[,i])]
    
#     # Create a dataframe for this method
#     method_data <- hallmarks %>%
#     filter(gene_symbol %in% input_genes) %>%  # Filter for input genes
#     group_by(gs_name) %>%
#     summarise(
#       total_genes_in_ontology = n(),
#       percentage_overlap = (total_genes_in_ontology / length(input_genes)) * 100
#     ) %>%
#     mutate(method = method) %>%  # Add method information
#     arrange(desc(percentage_overlap))
  
#     # Add this method's data to the results list
#     results[[method]] <- method_data
# }

# combined_stats <- bind_rows(results)

# # Plot the data with ontologies on the y-axis and methods on the x-axis, split into 3 panels
# ggplot(combined_stats, aes(x = method, y = gs_name, color = method, size = percentage_overlap)) +
#   geom_point(alpha = 0.7) +
#   scale_size_continuous(name = "Percentage of Overlap", range = c(2, 10)) + 
#   scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +  # Adjust spacing if needed
#   coord_cartesian(clip = "off") +  # Allow plot elements to extend beyond plot area
#   labs(
#     x = "Normalization Method",
#     y = "Hallmark",
#     title = "Percentage of detected HVGs in Hallmarks"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
#     plot.title = element_text(hjust = 0.5),             # Center the plot title
#     plot.margin = margin(t = 20, r = 20, b = 50, l = 20),  # Adjust margins for space
#     panel.background = element_rect(fill = "white", color = NA),  # White panel background
#     plot.background = element_rect(fill = "white", color = NA),   # White plot background
#     legend.background = element_rect(fill = "white", color = NA)  # White legend background
#   ) +
#   facet_wrap(~ gs_name, scales = "free_y", ncol = 1) +  # Split into 3 panels
#   theme(strip.text = element_text(size = 8))  # Adjust font size of panel labels

# # Save the plot and add filename
# ggsave(paste0(output_plot_path, '/Immune_analysis_facet.png'), width = 15, height = 50, limitsize = FALSE)
