args <- commandArgs(trailingOnly = TRUE)
color_list_path <- args[1]
cellTypeCSV <- args[2:length(args)]

# Load all cell types CSV
unique_cell_types <- c()
for (csv in cellTypeCSV) {
    cellType <- read.csv(csv)
    unique_cell_types <- unique(c(unique_cell_types, cellType$sup.clust))
}

unique_cell_types <- unique(unique_cell_types)
color_list <- read.csv(color_list_path, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
integers <- seq(from=1, to=(length(color_list[,1])-1), by=1)
random_colors <- sample(integers, length(unique_cell_types))
# random_colors <- seq(from=1, to=length(unique_cell_types), by=1)
# random_colors <- sample(c(0, 1, 2, 5, 11, 15, 20, 27, 29, 30, 32, 45, 49, 53, 74), length(unique_cell_types))

color_palette <- data.frame(cell_type = unique_cell_types, color = random_colors)

write.csv(color_palette, 'color_palette.csv', row.names = FALSE)