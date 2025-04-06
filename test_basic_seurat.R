suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})

work_dir = "projects/def-wainberg/karbabi/sc-benchmarking"
data_dir = "single-cell/SEAAD/subsampled"
source(file.path(work_dir, "utils_local.R"))

system_info()

size = "20K"

timers = TimerCollection(silent = FALSE)

# timers$with_timer("Load data (10X mtx)", {
#   data <- Read10X(
#       data.dir = paste0(data_dir, "/SEAAD_raw_", size),
#       gene.column=1)
#   data <- CreateSeuratObject(counts = data)
# })
# rm(data); invisible(gc())

# timers$with_timer("Load data (h5)", {
#   data = Read10X_h5(
#       filename = paste0(data_dir, "/SEAAD_raw_", size, ".h5"))
#   data <- CreateSeuratObject(counts = data)
# })
# rm(data); invisible(gc())

timers$with_timer("Load data", {
  data <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
})

# Note: QC filters are matched across libraries for timing, then 
# standardized by filtering to single_cell.py QC cells, not timed 

timers$with_timer("QC data", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
})

data <- subset(data, subset = passed_QC_tmp == TRUE)
print(paste0('cells: ', ncol(data), ', genes: ', nrow(data)))

# Note: No doublet detection

timers$with_timer("Normalize data", {
  data <- NormalizeData(data)
})

timers$with_timer("Find variable features", {
  data <- FindVariableFeatures(data)
})

timers$with_timer("PCA", {
  data = ScaleData(data)
  data = RunPCA(data)
})

timers$with_timer("Neighbor graph", {
  data = FindNeighbors(data)
})

timers$with_timer("Clustering (3 resolutions)", {
  for (resolution in c(0.5, 2, 1)) {
    data = FindClusters(data, resolution = resolution)
  }
})

print(paste0('seuratclusters: ', length(unique(data$seurat_clusters))))

# Note: There is no default number of PCs for UMAP

timers$with_timer("Embedding", {
  data = RunUMAP(data, dims = 1:10)
})

timers$with_timer("Plot embeddings", {
  DimPlot(data, reduction = "umap")
  # ggsave(paste0(work_dir, "/figures/seurat_embedding_cluster_", size, ".png"))
})

timers$with_timer("Find markers", {
  markers = FindAllMarkers(data)
})



timers$print_summary()
results_df = timers$to_dataframe()
print(results_df)





timers = TimerCollection(silent = FALSE)
timers$with_timer('load data', {
  df = data.frame(x = 1:10, y = rnorm(10))
})
timers$with_timer('calculate mean', {
  mean_val = mean(df$y)
})
timers$with_timer('filter values', {
  filtered = df[df$y > 0, ]
})
timers$print_summary()
results_df = timers$to_dataframe()
print(results_df)