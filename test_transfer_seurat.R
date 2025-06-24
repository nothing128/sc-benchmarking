suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
})

work_dir = "projects/sc-benchmarking"
data_dir = "~/single-cell/SEAAD"
source("utils_local.R")

args = commandArgs(trailingOnly=TRUE)
size <- args[1]
output <- args[2]

size_ref = c('1.2M' = '600K', '400K' = '200K', '20K' = '10K')

system_info()
timers = TimerMemoryCollection(silent = TRUE)

# Load data (query) ####
timers$with_timer("Load data (query)", {
  data_query <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
})

# Load data (ref) ####
timers$with_timer("Load data (ref)", {
  data_ref <- readRDS(paste0(data_dir, "/SEAAD_ref_", size_ref[size], ".rds"))
})

# Note: Temporaily add QC metrics 
# Running `prep_data.py` on the new cluster
# will add these metrics to a v5 seurat object

data_query[["nCount_RNA"]] <- colSums(data_query@assays$RNA@counts)
data_query[["nFeature_RNA"]] <- colSums(data_query@assays$RNA@counts > 0)

# Quality control ####
timers$with_timer("Quality control", {
  data_query[["percent.mt"]] <- PercentageFeatureSet(data_query, pattern = "^MT-")
  data_query <- subset(data_query, subset = nFeature_RNA > 200 & percent.mt < 5)
})

# doublet detection (using doubletFinder)
timers$with_timer("Doublet detection", {
  # preprocessing for Doublet detection 
  data_query <- NormalizeData(data_query)
  data <- FindVariableFeatures(data_query)
  data_query <- ScaleData(data_query)
  data_query <- RunPCA(data_query)
  data_query <- RunUMAP(data_query, dims = 1:10)
  # pK identification (no ground truth)
  sweep.res.data_query <- paramSweep(data_query, PCs = 1:10, sct = FALSE)
  sweep.stats_dat_querya <- summarizeSweep(sweep.res.data_query, GT = FALSE)
  bcmvn <- find.pK(sweep.stats_data_query)
  homotypic.prop <- modelHomotypic(data_query@meta.data$seurat_clusters)
  # Homotypic Doublet Proportion Estimate
  n_cells <- ncol(data_query)
  doublet_rate <- (n_cells / 1000) * 0.008          
  nExp_poi <- round(doublet_rate* n_cells)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # detection
  data_query <- doubletFinder(data_query, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

})

# remove doublets (not timed)
df_classification_col <- names(data_query@meta.data)[grepl("DF.classifications", names(data_query@meta.data))]
data_query <- subset(data_query, subset = !!sym(df_classification_col) == "Singlet")

# Normalization ####
timers$with_timer("Normalization", {
  data_ref <- NormalizeData(data_ref)
  data_query <- NormalizeData(data_query)
})

# Feature selection ####
timers$with_timer("Feature selection", {
  data_ref <- FindVariableFeatures(data_ref)
  data_query <- FindVariableFeatures(data_query)
})

# PCA ####
timers$with_timer("PCA", {
  data_ref <- ScaleData(data_ref)
  data_ref <- RunPCA(data_ref)
  data_query <- ScaleData(data_query)
  data_query <- RunPCA(data_query)
})

# Transfer labels ####
timers$with_timer("Transfer labels", {
  anchors <- FindTransferAnchors(
    reference = data_ref, query = data_query, dims = 1:30,
    reference.reduction = "pca")
  predictions <- TransferData(
    anchorset = anchors, refdata = data_ref$subclass)
  data_query <- AddMetaData(
    object = data_query, metadata = predictions)
})

print("--- Transfer Accuracy ---")
data_query$prediction.match <- data_query$predicted.id == data_query$subclass
result <- data_query@meta.data %>%
  group_by(subclass) %>%
  summarize(
    n_correct = sum(prediction.match),
    n_total = n()
  ) %>%
  ungroup() %>%
  bind_rows(
    summarize(
      .,
      subclass = "Total",
      n_correct = sum(n_correct),
      n_total = sum(n_total)
    )
  ) %>%
  mutate(percent_correct = (n_correct / n_total) * 100) 

accuracy_path <- paste0(work_dir,"/output/test_transfer_seurat_", size,"_accuracy.csv")
write.csv(result, accuracy_path, row.names = FALSE)
print(result,n = Inf)

timers$print_summary(sort = FALSE)

timers_df <- timers$to_dataframe(unit = "s", sort = FALSE) %>%
  mutate(
    test = 'test_transfer_seurat',
    size = size
  )

write.csv(timers_df, output, row.names = FALSE)
