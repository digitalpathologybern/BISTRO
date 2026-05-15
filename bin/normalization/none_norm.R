suppressPackageStartupMessages({
  library(data.table)
})


# Read terminal parameters
args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
output_csv_path <- args[2]

# Load count matrix and metadata csv
# counts <- read.csv(count_matrix_path, row.names = 1)
message("Loading count matrix...")
dt <- fread(count_matrix_path, check.names = FALSE)
cellnames <- dt[[1]]
dt[[1]] <- NULL
counts <- as.matrix(dt)
rownames(counts) <- as.character(cellnames)
rm(dt); gc()
message("Counts loaded (cells x genes): ", paste(dim(counts), collapse = " x "))

# Get a vector of 1 for each column
sf <- rep(1, nrow(counts))
# Get input name and add _cpm 
filename <- strsplit(count_matrix_path, split = "/")[[1]][length(strsplit(count_matrix_path, split = "/")[[1]])]
filename <- strsplit(toString(filename), split='.csv')

df <- data.frame(sf)
colnames(df) <- c("x")
rownames(df) <- rownames(counts)

write.table(counts, 
            paste0(output_csv_path, '/', filename, '_none.csv'), 
            sep = ',', row.names = T, col.names = NA, quote = F)

write.table(df, 
            paste0(output_csv_path, '/', filename, '_noneSF.csv'), 
            sep = ',', row.names = T, col.names = NA, quote = F)
