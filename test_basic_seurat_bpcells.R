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

system_info()
timers = TimerCollection(silent = TRUE)

# Note: Loading is much slower from $SCRATCH disk

# Not timed
if (file.exists(paste0(data_dir, "/bpcells/", size))) {
  unlink(paste0(data_dir, "/bpcells/", size), recursive = TRUE)
}

timers$with_timer("Load data", {
    mat_disk <- open_matrix_10x_hdf5(
      path = paste0(data_dir, "/SEAAD_raw_", size,".h5"))
    mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")

    file_path = paste0(data_dir, "/bpcells/", size)
    write_matrix_dir(
      mat = mat_disk,
      dir = file_path
    )
    mat <- open_matrix_dir(dir = file_path)
    data <- CreateSeuratObject(counts = mat)
  })

# Not timed
data <- AddMetaData(data, metadata = data.frame(
  nCount_RNA = colSums(data@assays$RNA@layers$counts),
  nFeature_RNA = colSums(data@assays$RNA@layers$counts > 0)
))

timers$with_timer("Quality control", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
})

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

# svd <- BPCells::svds(data@assays$RNA@layers$data, k=50)
# pca <- multiply_cols(svd$v, svd$d)
# data[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_")

timers$with_timer("Neighbor graph", {
  data = FindNeighbors(data, dims = 1:10)
})

#TODO: The number of clusters needs to match across libraries

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
timers_df$test = 'test_basic_seurat_bpcells'
timers_df$size = size

write.csv(timers_df, output, row.names = FALSE)
unlink(paste0(data_dir, "/bpcells/", size), recursive = TRUE)
