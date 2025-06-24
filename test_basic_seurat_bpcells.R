suppressPackageStartupMessages({
    library(tidyverse)
    library(BPCells)
    library(Seurat)
    library(DoubletFinder)
})

work_dir = "projects/sc-benchmarking"
data_dir = "~/single-cell/SEAAD"
source("utils_local.R")

args = commandArgs(trailingOnly=TRUE)
size <- args[1]
output <- args[2]

scratch_dir <- Sys.getenv("SCRATCH")
bpcells_dir <- file.path(scratch_dir, "bpcells")
if (!dir.exists(bpcells_dir)) {
    dir.create(bpcells_dir, recursive = TRUE)
}

system_info()
timers = TimerMemoryCollection(silent = TRUE)

# Not timed
if (file.exists(file.path(bpcells_dir, size))) {
  unlink(file.path(bpcells_dir, size), recursive = TRUE)
}

# Load data ####
timers$with_timer("Load data", {
    mat_disk <- open_matrix_10x_hdf5(
      path = paste0(data_dir, "/SEAAD_raw_", size,".h5"))
    mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")

    file_path = paste0(bpcells_dir, "/bpcells/", size)
    write_matrix_dir(
      mat = mat_disk,
      dir = file_path
    )
    mat <- open_matrix_dir(dir = file_path)
    data <- CreateSeuratObject(counts = mat)
  })

# Not timed
# data <- AddMetaData(data, metadata = data.frame(
#   nCount_RNA = colSums(data@assays$RNA@layers$counts),
#   nFeature_RNA = colSums(data@assays$RNA@layers$counts > 0)
# ))
data[["nCount_RNA"]] <- colSums(data@assays$RNA@counts)
data[["nFeature_RNA"]] <- colSums(data@assays$RNA@counts > 0)

# Quality control ####
timers$with_timer("Quality control", {
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
})

# doublet detection (using doubletFinder)
timers$with_timer("Doublet detection", {
  # preprocessing for Doublet detection 
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data)
  data <- ScaleData(data)
  data <- RunPCA(data)
  data <- RunUMAP(data, dims = 1:10)
  # pK identification (no ground truth)
  sweep.res.data <- paramSweep(data, PCs = 1:10, sct = FALSE)
  sweep.stats_data <- summarizeSweep(sweep.res.data, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_data)
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  # Homotypic Doublet Proportion Estimate
  n_cells <- ncol(data)
  doublet_rate <- (n_cells / 1000) * 0.008          
  nExp_poi <- round(doublet_rate* n_cells)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # detection
  data <- doubletFinder(data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

})

# remove doublets (not timed)
df_classification_col <- names(data@meta.data)[grepl("DF.classifications", names(data@meta.data))]
data <- subset(data, subset = !!sym(df_classification_col) == "Singlet")

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

# svd <- BPCells::svds(data@assays$RNA@layers$data, k=50)
# pca <- multiply_cols(svd$v, svd$d)
# data[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_")

# Neighbor graph ####
timers$with_timer("Neighbor graph", {
  data = FindNeighbors(data, dims = 1:10)
})

#TODO: The number of clusters needs to match across libraries

# Clustering (3 resolutions) ####
timers$with_timer("Clustering (3 resolutions)", {
  for (resolution in c(0.5, 2, 1)) {
    data = FindClusters(data, resolution = resolution)
  }
})

print(paste0('seuratclusters: ', length(unique(data$seurat_clusters))))

# Embedding ####
timers$with_timer("Embedding", {
  data = RunUMAP(data, dims = 1:10)
})

# Plot embeddings ####
timers$with_timer("Plot embeddings", {
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
timers_df$test = 'test_basic_seurat_bpcells'
timers_df$size = size

write.csv(timers_df, output, row.names = FALSE)
unlink(paste0(bpcells_dir, "/bpcells/", size), recursive = TRUE)
