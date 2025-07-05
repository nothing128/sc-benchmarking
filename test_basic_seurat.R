suppressPackageStartupMessages({
  library(tidyverse)
  library(BPCells)
  library(Seurat)
})

work_dir = "projects/sc-benchmarking"
data_dir = "single-cell/SEAAD"
source(file.path(work_dir, "utils_local.R"))

args = commandArgs(trailingOnly=TRUE)
size <- args[1]
output <- args[2]

scratch_dir <- Sys.getenv("SCRATCH")
bpcells_dir_test <- file.path(scratch_dir, "bpcells", "basic")
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
    path = paste0(data_dir, "/SEAAD_raw_", size,".h5ad"))
  mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")
  file_path = file.path(bpcells_dir_test, size)
  write_matrix_dir(
    mat = mat_disk,
    dir = file_path
  )
  mat <- open_matrix_dir(dir = file_path)
  # Custom utility to read obs metadata from h5ad file
  obs_metadata <- read_h5ad_obs(paste0(data_dir, "/SEAAD_raw_", size,".h5ad"))
  data <- CreateSeuratObject(counts = mat, meta.data = obs_metadata)
})

# Quality control ####
timers$with_timer("Quality control", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
})

# Note: No doublet detection offered in Seurat

# Normalization ####
timers$with_timer("Normalization", {
  data <- NormalizeData(
    data, normalization.method = "LogNormalize", scale.factor = 10000)
})

# Feature selection ####
timers$with_timer("Feature selection", {
  data <- FindVariableFeatures(
    data, selection.method = "vst", nfeatures = 2000)  
})

# PCA ####
timers$with_timer("PCA", {
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  data <- RunPCA(data, features = VariableFeatures(object = data))
})

# Neighbor graph ####
timers$with_timer("Nearest neighbors", {
  data = FindNeighbors(data, dims = 1:10)
})

# Clustering (3 res.) ####
timers$with_timer("Clustering (3 res.)", {
  for (resolution in c(0.5, 2, 1)) {
    data = FindClusters(data, resolution = resolution)
  }
})

# Embedding ####
timers$with_timer("Embedding", {
  data = RunUMAP(data, dims = 1:10)
})

# Plot embeddings ####
timers$with_timer("Plot embedding", {
  DimPlot(data, reduction = "umap", group.by = "subclass")
  ggsave(paste0(work_dir, "/figures/seurat_embedding_subclass_", size, ".png"),
        dpi = 300, units = "in", width = 10, height = 10)
})

# Find markers ####
timers$with_timer("Find markers", {
  markers = FindAllMarkers(data, group.by = "subclass", only.pos = TRUE)
})

timers$print_summary(sort = FALSE)
timers_df = timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$library = 'seurat'
timers_df$test = 'basic'
timers_df$size = size

write.csv(timers_df, output, row.names = FALSE)

unlink(file.path(bpcells_dir_test, size), recursive = TRUE)
rm(data, markers, timers, timers_df, mat, mat_disk, obs_metadata)
gc()
