library(SpaNorm)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(Matrix)
library(data.table)
library(BiocParallel)

set.seed(0)


# ============================================================================
# NOTE: FOV rasterization has been removed from this script.
# FOV assignment is now centralized in bin/preprocessing/assign_fov.py
# and runs as a Nextflow process (assignFOV) before normalization.
# The 'fov' column is already present in the enriched metadata CSV.
#
# If you need to modify rasterization logic, edit:
#   utils/helpers.py -> rasterize_fov_grid()
#   bin/preprocessing/assign_fov.py
# ============================================================================


# Read terminal parameters
# Args: count_matrix, metadata, output_dir, separate_fovs, technology, transformation
# NOTE: 'rasterize' parameter has been removed (was args[6]).
#       Transformation is now args[6] instead of args[7].
args <- commandArgs(trailingOnly = TRUE)
count_matrix_path <- args[1]
metadata_matrix_path <- args[2]
output_csv_path <- args[3]
separate_fovs <- args[4]
technology <- args[5]
transformation <- args[6]

# Load count matrix
dt <- fread(count_matrix_path, check.names = FALSE)
cellnames <- dt[[1]]
dt[[1]] <- NULL

# Convert to numeric in-place (no copy)
for (j in seq_len(ncol(dt))) set(dt, j = j, value = as.numeric(dt[[j]]))

# Create a sparse matrix directly without transposing yet
counts <- as(as.matrix(dt), "dgCMatrix")
rownames(counts) <- as.character(cellnames)
rm(dt); gc()

# Transpose lazily (avoid double memory)
counts <- Matrix::t(counts)
gc()
message("Counts loaded: ", nrow(counts), " genes x ", ncol(counts), " cells")

# Load metadata (FOV column is already present from assign_fov.py)
meta <- read.csv(metadata_matrix_path)
meta_rownames <- read.csv(
  metadata_matrix_path, 
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = c("character")
)
rownames(meta) <- meta_rownames[,1]
meta <- meta[,-1, drop = FALSE]

# Verify FOV column exists
if (!"fov" %in% colnames(meta)) {
  stop("ERROR: 'fov' column not found in metadata. Ensure assignFOV process has run.")
}

message("FOV identifiers loaded from metadata: ", length(unique(meta$fov)), " unique FOVs")


# Get input name 
filename <- strsplit(count_matrix_path, split = "/")[[1]][length(strsplit(count_matrix_path, split = "/")[[1]])]
filename <- strsplit(toString(filename), split='.csv')

# Create SpatialExperiment object
# All platforms now use x_global_um / y_global_um or x_local_um / y_local_um
# (standardized um-unit coordinate schema)
if (technology == 'CosMx') {
    spe <- SpatialExperiment(
        assay = list(counts = counts),
        colData = meta,
        spatialCoordsNames = c("x_global_um", "y_global_um")
    )
    sampling <- 0.3
} else if (technology == 'Xenium') {
    spe <- SpatialExperiment(
        assay = list(counts = counts),
        colData = meta,
        spatialCoordsNames = c("x_local_um", "y_local_um")
    )
    sampling <- 1
} else if (technology == 'MERFISH') {
    spe <- SpatialExperiment(
        assay = list(counts = counts),
        colData = meta,
        spatialCoordsNames = c("x_local_um", "y_local_um")
    )
    sampling <- 1
}

# Precompute size factors
spe <- fastSizeFactors(spe)

# ============================================================================
# SpaNorm-specific FOV merging: merge small FOVs (<300 cells) into the
# preceding FOV to avoid singularities in spline fitting.
# This merging is NOT applied to the pipeline-wide FOV assignments
# (handled by assign_fov.py with no merging), only to SpaNorm's internal
# per-FOV processing.
# ============================================================================
if (separate_fovs == "1") {
    min_cells_spanorm <- 300
    fov_counts <- as.data.frame(table(meta$fov))
    colnames(fov_counts) <- c("fov_id", "cell_count")
    fov_counts$fov_id <- as.integer(as.character(fov_counts$fov_id))
    fov_counts <- fov_counts[order(fov_counts$fov_id), ]

    fov_remap <- list()
    new_fov <- 0
    prev_valid <- 0

    for (i in seq_len(nrow(fov_counts))) {
        fov_id <- fov_counts$fov_id[i]
        count <- fov_counts$cell_count[i]
        if (i == 1 || count >= min_cells_spanorm) {
            fov_remap[[as.character(fov_id)]] <- new_fov
            prev_valid <- new_fov
            new_fov <- new_fov + 1
        } else {
            fov_remap[[as.character(fov_id)]] <- prev_valid
        }
    }

    original_n_fovs <- length(unique(spe$fov))
    spe$fov <- sapply(as.character(spe$fov), function(x) fov_remap[[x]])
    merged_n_fovs <- length(unique(spe$fov))
    message("SpaNorm FOV merging: ", original_n_fovs, " -> ", merged_n_fovs,
            " FOVs (min_cells=", min_cells_spanorm, ")")

    assay(spe, 'logcounts') <- counts(spe)

    # Get unique FOVs (from the merged fov column)
    fov_list <- sort(unique(spe$fov))
    message("Processing ", length(fov_list), " FOVs separately")

    # 1. Prepare: Get the full list of column indices sorted by FOV
    all_cell_indices <- seq_len(ncol(spe))
    cell_order_map <- order(spe$fov) 
    original_col_indices <- all_cell_indices[cell_order_map]
    combined_logcounts_matrix <- matrix(NA, nrow=nrow(spe), ncol=ncol(spe)) 

    ############################################################################################################
    # CHECKPOINT VERSION
    ############################################################################################################
    # 2. Process FOVs with per-FOV checkpointing
    pub_dir <- Sys.getenv("SPANORM_CHECKPOINT_DIR", output_csv_path)
    ckpt_dir <- file.path(pub_dir, "spanorm_checkpoints")
    dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

    logcounts_list <- list()

    for (i in seq_along(fov_list)) {
        fov <- fov_list[i]
        ckpt_file <- file.path(ckpt_dir, paste0("fov_", fov, ".rds"))

        # Skip if checkpoint exists
        if (file.exists(ckpt_file)) {
            cat(sprintf("[%d/%d] FOV %s: loading from checkpoint\n",
                        i, length(fov_list), fov))
            flush.console()
            logcounts_list[[as.character(fov)]] <- readRDS(ckpt_file)
            next
        }

        cat(sprintf("[%d/%d] FOV %s: processing (%d cells) ...\n",
                    i, length(fov_list), fov, sum(spe$fov == fov)))
        flush.console()

        fov_cols <- spe$fov == fov
        current_counts <- spe[, fov_cols]

        if (sum(fov_cols) < 1000) {
            current_counts <- SpaNorm(current_counts, sample.p = sampling,
                                      df.tps = 3, adj.method = transformation)
        } else {
            current_counts <- tryCatch({
                SpaNorm(current_counts, sample.p = sampling,
                        adj.method = transformation)
            }, error = function(e) {
                cat(sprintf("  Retrying FOV %s with df.tps = 3\n", fov))
                flush.console()
                SpaNorm(current_counts, sample.p = sampling,
                        df.tps = 3, adj.method = transformation)
            })
        }

        normalized_logcounts <- assay(current_counts, "logcounts")

        # Save checkpoint
        saveRDS(normalized_logcounts, ckpt_file)
        cat(sprintf("  Checkpoint saved: %s\n", basename(ckpt_file)))
        flush.console()

        logcounts_list[[as.character(fov)]] <- normalized_logcounts
        rm(current_counts, normalized_logcounts); gc()
    }

    # Clean up checkpoint dir after all FOVs complete
    unlink(ckpt_dir, recursive = TRUE)
    cat("All FOVs complete, checkpoints cleaned up\n")
    flush.console()



    ############################################################################################################
    # ORIGINAL VERSION
    ############################################################################################################
    # # 2. Initialize a list for the normalized results
    # logcounts_list <- list() 

    # for(fov in fov_list) {
    #   message('================================================')
    #   message(fov)
    #   message('================================================')

    #   # Identify the columns for the current FOV
    #   fov_cols <- spe$fov == fov 
          
    #   current_counts <- spe[, fov_cols]

    #   # If we have less than 1000 cells, reduce splines df to 3
    #   if (sum(fov_cols) < 1000) {
    #     current_counts <- SpaNorm(current_counts, sample.p = sampling, df.tps = 3, adj.method = transformation) 
    #   } else {
    #     current_counts <- tryCatch({
    #       SpaNorm(current_counts, sample.p = sampling)
    #     }, error = function(e) {
    #       print(paste("Error encountered with default df.tps for FOV", fov, ". Retrying with df.tps = 3."))
    #       return(SpaNorm(current_counts, sample.p = sampling, df.tps = 3, adj.method = transformation))
    #     })
    #   }

    #   normalized_logcounts <- assay(current_counts, 'logcounts')
    #   print(normalized_logcounts[1:5, 1:5])
          
    #   # 3. Store the logcounts matrix
    #   logcounts_list[[as.character(fov)]] <- normalized_logcounts
    # }

    # # # # n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "6"))
    # # # # logcounts_list <- bplapply(fov_list, function(fov) {
    # # # #     message('================================================')
    # # # #     message(fov)
    # # # #     message('================================================')
    # # # #     fov_cols <- spe$fov == fov
    # # # #     current_spe <- spe[, fov_cols]

    # # # #     if (sum(fov_cols) < 1000) {
    # # # #         current_spe <- SpaNorm(current_spe, sample.p = sampling,
    # # # #                               df.tps = 3, adj.method = transformation)
    # # # #     } else {
    # # # #         current_spe <- tryCatch({
    # # # #             SpaNorm(current_spe, sample.p = sampling, adj.method = transformation)
    # # # #         }, error = function(e) {
    # # # #             message("Retrying FOV ", fov, " with df.tps = 3")
    # # # #             SpaNorm(current_spe, sample.p = sampling,
    # # # #                     df.tps = 3, adj.method = transformation)
    # # # #         })
    # # # #     }

    # # # #     assay(current_spe, "logcounts")
    # # # # }, BPPARAM = MulticoreParam(workers = n_cores))
    # # # # names(logcounts_list) <- as.character(fov_list)

    # 4. Final Assignment: Merge and Write ONCE
    combined_matrix <- do.call(cbind, logcounts_list[as.character(unique(spe$fov))])

    original_colnames <- colnames(spe)
    combined_colnames <- colnames(combined_matrix)

    match_indices <- match(original_colnames, combined_colnames)
    final_logcounts_matrix <- combined_matrix[, match_indices]

    # Assign to the main object
    assay(spe, 'logcounts') <- final_logcounts_matrix 
    print("Completed assignment of all logcounts to spe object.")
    
} else {

spe <- SpaNorm(spe, sample.p = 0.35, adj.method = transformation)

}

# Save data matrix
write.table(t(as.matrix(logcounts(spe))), 
            paste0(output_csv_path, '/', filename, '_spanorm-', transformation ,'.csv'),
            sep = ',', 
            row.names = TRUE,
            col.names = NA,
            quote = FALSE)