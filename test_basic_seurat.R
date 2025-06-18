suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(Seurat)
})  

work_dir = "projects/sc-benchmarking"
data_dir = "~/single-cell/SEAAD"
source(file.path("utils_local.R"))

args = commandArgs(trailingOnly=TRUE)
size <- args[1]
output <- args[2]

system_info()
timers = TimerCollection(silent = TRUE)

# Note: Loading is much slower from $SCRATCH disk

timers$with_timer("Load data", {
  data <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
})

# Note: Temporaily add QC metrics 
# Running `prep_data.py` on the new cluster
# will add these metrics to a v5 seurat object

data[["nCount_RNA"]] <- colSums(data@assays$RNA@counts)
data[["nFeature_RNA"]] <- colSums(data@assays$RNA@counts > 0)

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

timers$with_timer("Neighbor graph", {
  data <- FindNeighbors(data, dims = 1:10)
})

#TODO: The number of clusters needs to match across libraries

timers$with_timer("Clustering (3 resolutions)", {
  for (resolution in c(0.5, 2, 1)) {
    data <- FindClusters(data, resolution = resolution)
  }
})

print(paste0('seuratclusters: ', length(unique(data$seurat_clusters))))

timers$with_timer("Embedding", {
  data <- RunUMAP(data, dims = 1:10)
})

timers$with_timer("Plot embeddings", {
  DimPlot(data, reduction = "umap", group.by = "subclass")
  ggsave(paste0(work_dir, "/figures/seurat_embedding_subclass_", size, ".png"),
        dpi = 300, units = "in", width = 10, height = 10)
})

timers$with_timer("Find markers", {
  markers <- FindAllMarkers(data, group.by = "subclass", only.pos = TRUE)
})

timers$print_summary(sort = FALSE)
timers_df <- timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$test <- 'test_basic_seurat'
timers_df$size <- size

write.csv(timers_df, output, row.names = FALSE)
