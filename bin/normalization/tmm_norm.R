suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
  library(data.table)
})

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
counts <- t(counts)  # Transpose to genes x cells
rm(dt); gc()
message("Counts loaded (genes x cells): ", paste(dim(counts), collapse = " x "))

# Get input name and add _cpm 
filename <- strsplit(count_matrix_path, split = "/")[[1]][length(strsplit(count_matrix_path, split = "/")[[1]])]
filename <- strsplit(toString(filename), split='.csv')

# Get the scaling values and calculate the effective library size
scaling_values <- normLibSizes(counts, method='TMM')
print('Scaling values calculated')
effective_lib_size <- colSums(counts) * scaling_values

# Apply the effective library size to the count data
TMM <- sweep(counts, 2, effective_lib_size, FUN = "/")

# Save table with the correct format
tmm_table <- data.frame(t(TMM))
write.table(tmm_table, 
            paste0(output_csv_path, '/', filename, '_tmm.csv'), 
            sep = ',', row.names = T, col.names = NA, quote = F)
print('TMM normalized counts saved')
rm(TMM); gc()

# Transform effective_lib_size to data frame with rownames as sample names and colnames as x
df <- data.frame(effective_lib_size)
colnames(df) <- 'x'
rownames(df) <- colnames(counts)
write.table(df, 
            paste0(output_csv_path, '/', filename, '_tmmSF.csv'), 
            sep = ',', row.names = T, col.names = NA, quote = F)