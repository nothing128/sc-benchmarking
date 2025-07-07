suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(Seurat)
  library(BPCells)
})  

work_dir = "projects/sc-benchmarking"
data_dir = "~/single-cell/SEAAD"
source(file.path("utils_local.R"))

args = commandArgs(trailingOnly=TRUE)
size <- args[1]
output <- args[2]

scratch_dir <- Sys.getenv("SCRATCH")
bpcells_dir_test <- file.path(scratch_dir, "bpcells", "de")
if (!dir.exists(bpcells_dir_test)) {
    dir.create(bpcells_dir_test, recursive = TRUE)
}

system_info()
timers = TimerMemoryCollection(silent = TRUE)

# Not timed
if (file.exists(file.path(bpcells_dir_test, size))) {
  unlink(file.path(bpcells_dir_test, size), recursive = TRUE)
}

# Load data ####
timers$with_timer("Load data", {
    mat_disk <- open_matrix_anndata_hdf5(
      path = file.path(data_dir, paste0("SEAAD_raw_", size,".h5ad")))
    mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")
    file_path <- file.path(bpcells_dir_test, size)
    write_matrix_dir(
      mat = mat_disk,
      dir = file_path
    )
    mat <- open_matrix_dir(dir = file_path)
    # Custom utility to read obs metadata from h5ad file
    obs_metadata <- read_h5ad_obs(
      file.path(data_dir, paste0("SEAAD_raw_", size,".h5ad")))
    data <- CreateSeuratObject(counts = mat, meta.data = obs_metadata)
  })


# Quality control ####
timers$with_timer("Quality control", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
})

# TODO: Add doublet detection
# Subset the BPCells matrix to each sample and loop through each sample

# Not timed
data$ad_dx <- ifelse(data$ad_dx == "1", "AD", "Control")

# Data transformation (pseudobulk / normalization) ####
timers$with_timer("Data transformation (pseudobulk / normalization)", {
  data <- AggregateExpression(
    data, 
    assays = "RNA", 
    return.seurat = TRUE, 
    group.by = c("ad_dx", "sample", "subclass")
  )
})

# Differential expression ####
timers$with_timer("Differential expression", {
  data$group <- paste(data$subclass, data$ad_dx, sep = "_")
  Idents(data) <- "group"

  de_list <- list()
  for (subclass in unique(data$subclass)) {
    de_list[[subclass]] <- FindMarkers(
      data, 
      ident.1 = paste(subclass, "AD", sep = "_"), 
      ident.2 = paste(subclass, "Control", sep = "_"), 
      test.use = "DESeq2"
    )
  }
  de <- do.call(rbind, de_list)
})

timers$print_summary(sort = FALSE)
<<<<<<< HEAD
timers_df <- timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$test <- 'test_de_seurat'
timers_df$size <- size

=======
>>>>>>> 477e7cfb0cb6f0253d124874bf13f7e15b546892
timers_df <- timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$library <- 'seurat'
timers_df$test <- 'de'
timers_df$size <- size

write.csv(timers_df, output, row.names = FALSE)

unlink(file.path(bpcells_dir_test, size), recursive = TRUE)
rm(data, de, de_list, timers, timers_df, mat, mat_disk, obs_metadata)
gc()


# # Normalization ####
# timers$with_timer("Normalization", {
#   data <- NormalizeData(
#     data, normalization.method = "LogNormalize", scale.factor = 10000)
# })

# # Differential expression ####
# timers$with_timer("Differential expression (wilcoxon)", {
#   data$group <- paste(data$subclass, data$ad_dx, sep = "_")
#   Idents(data) <- "group"

#   de_list <- list()
#   for (subclass in unique(data$subclass)) {
#     de_list[[subclass]] <- FindMarkers(
#       data, 
#       ident.1 = paste(subclass, "AD", sep = "_"), 
#       ident.2 = paste(subclass, "Control", sep = "_")
#     )
#   }
#   de <- do.call(rbind, de_list)
# })