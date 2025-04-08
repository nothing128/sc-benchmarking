suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})

work_dir = "projects/sc-benchmarking"
data_dir = "scratch/single-cell/SEAAD"
source(file.path(work_dir, "utils_local.R"))

system_info()

size = "20K"

timers = TimerCollection(silent = FALSE)

# # Note: Load times should not be considered when using $SCRATCH disk 

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

# # Note: rds files contains additional metadata columns vs above

timers$with_timer("Load data", {
  data <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
})

# Note: QC filters are matched across libraries for timing, then 
# standardized by filtering to single_cell.py QC cells, not timed 

timers$with_timer("Quality control", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
})

data <- subset(data, subset = passed_QC_tmp == TRUE)
print(paste0('cells: ', ncol(data), ', genes: ', nrow(data)))

# Note: No doublet detection offered in Seurat

timers$with_timer("Normalization", {
  data <- NormalizeData(
    data, normalization.method = "LogNormalize", scale.factor = 10000)
})

timers$with_timer("Feature selection", {
  data <- FindVariableFeatures(
    data, selection.method = "vst", nfeatures = 2000)  
})

timers$with_timer("PCA", {
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  data <- RunPCA(data, features = VariableFeatures(object = data))
})

timers$with_timer("Neighbor graph", {
  data = FindNeighbors(data, dims = 1:10)
})

timers$with_timer("Clustering (3 resolutions)", {
  for (resolution in c(0.5, 2, 1)) {
    data = FindClusters(data, resolution = resolution)
  }
})

print(paste0('seuratclusters: ', length(unique(data$seurat_clusters))))

timers$with_timer("Embedding", {
  data = RunUMAP(data, dims = 1:10)
})

timers$with_timer("Plot embeddings", {
  DimPlot(data, reduction = "umap")
  ggsave(paste0(work_dir, "/figures/seurat_embedding_cluster_", size, ".png"),
         dpi = 300, units = "in", width = 10, height = 10)
})

timers$with_timer("Find markers", {
  markers = FindAllMarkers(data, only.pos = TRUE)
})

timers$print_summary(sort = FALSE)
timers_df = timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$test = 'test_basic_seurat'
timers_df$size = size

print(timers_df)

write.csv(timers_df, 
  paste0(work_dir, "/output/test_basic_seurat_", size, ".csv"), 
  row.names = FALSE)
