suppressPackageStartupMessages({
    library(dplyr)
    library(BPCells)
    library(Seurat)
    library(SeuratObject)
    library(ggplot2)
    library(patchwork)
    library(hdf5r)
    library(Azimuth)
})

work_dir = "projects/rrg-wainberg/lamming6/sc-benchmarking"
data_dir = "single-cell/SEAAD/subsampled"
source(file.path("utils_local.R"))

system_info()

for (size in c("20K")) {  
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

  # timers$with_timer("Load data (h5Seurat)", {
  #   data <- LoadH5Seurat(
  #     data.dir = paste0(data_dir, "/SEAAD_raw_", size, ".h5Seurat"))
  # })
  # rm(data); invisible(gc())

  # Note: Loading is much slower from $SCRATCH disk

  timers$with_timer("Load data (h5ad/rds)", {
    mat_raw <- open_matrix_10x_hdf5(
      path = paste0(data_dir, "/SEAAD_raw_", size,".h5")
    )
    # Write the matrix to a directory
    write_matrix_dir(
      mat = mat_raw,
      dir = work_dir
    )
    # Now that we have the matrix on disk, we can load it
    mat <- open_matrix_dir(dir = work_dir)
    mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
    
    data <- CreateSeuratObject(counts = mat_raw)
  })

  # Add nFeautre_RNA to DS TODO: should this be part of QC
  # if ("nFeature_RNA" %in% colnames(data@meta.data)) {
  # print("Warning: nFeature_RNA already exists. It will be overwritten.")
  # }
  # assay_name <- "RNA" 
  # if (!assay_name %in% Assays(data)) {
  #   stop(paste("Assay '", assay_name, "' not found in the Seurat object. Available assays are: ",
  #             paste(Assays(data), collapse=", ")))
  # }
  # counts_matrix <- GetAssayData(object = data, assay = assay_name, layer = "counts")
  # nFeature_RNA_values <- Matrix::colSums(counts_matrix > 0)
  # data$nFeature_RNA <- nFeature_RNA_values


  # Note: QC filters are matched across libraries for timing, then 
  # standardized by filtering to single_cell.py QC cells, not timed 

  timers$with_timer("Quality control", {
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    print(summary(data$percent.mt)) # TODO: all NAs???
    data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 5)
    # data <- subset(data, subset = nFeature_RNA > 200)
  })

  # data <- subset(data, subset = tmp_passed_QC == TRUE)
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
    # num_features <- 50
    # svd <- BPCells::svds(data, k=num_features)
    # # Alternate option: irlba::irlba(mat_norm, nv=50)
    # data <- multiply_cols(svd$v, svd$d)
    data <- RunPCA(data, features = VariableFeatures(object = data))
  })

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

  timers$with_timer("Save data", {
    saveRDS(data, paste0(work_dir, "/output/test_write.rds"))
  })

  timers$print_summary(sort = FALSE)
  timers_df = timers$to_dataframe(unit = "s", sort = FALSE)
  timers_df$test = 'test_basic_seurat'
  timers_df$size = size

  print(timers_df)

  write.csv(timers_df, 
    paste0(work_dir, "/output/test_basic_seurat_", size, ".csv"), 
    row.names = FALSE)
}

'''
--- Timing Summary 20K ---
Load data (h5ad/rds) took 12s 254ms (4.3%)
Quality control took 6s 109ms (2.1%)
Normalization took 5s 102ms (1.8%)
Feature selection took 16s 516ms (5.8%)
PCA took 35s 672ms (12.5%)
Neighbor graph took 4s 817ms (1.7%)
Clustering (3 resolutions) took 4s 757ms (1.7%)
Embedding took 30s 45ms (10.5%)
Plot embeddings took 2s 801ms (1.0%)
Find markers took 2m 47s (58.6%)

Total time: 4m 45s

'''