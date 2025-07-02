suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(BPCells)
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

size_ref = c('1.2M' = '600K', '400K' = '200K', '20K' = '10K')

system_info()
timers = TimerMemoryCollection(silent = TRUE)

# Not timed
if (file.exists(file.path(bpcells_dir, size))) {
  unlink(file.path(bpcells_dir, size), recursive = TRUE)
}
if (file.exists(file.path(bpcells_dir, paste0("ref_", size_ref[size])))) {
  unlink(file.path(bpcells_dir, paste0("ref_", size_ref[size])), recursive = TRUE)
}

# Load data (query) ####
timers$with_timer("Load data (query)", {
  mat_disk <- open_matrix_anndata_hdf5(
    path = paste0(data_dir, "/SEAAD_raw_", size,".h5ad"))
  mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")
  file_path = file.path(bpcells_dir, size)
  write_matrix_dir(
    mat = mat_disk,
    dir = file_path
  )
  mat <- open_matrix_dir(dir = file_path)
  data_query <- CreateSeuratObject(counts = mat)
})

# Load data (ref) ####
timers$with_timer("Load data (ref)", {
  mat_disk <- open_matrix_anndata_hdf5(
    path = paste0(data_dir, "/SEAAD_ref_", size_ref[size],".h5ad"))
  mat_disk <- convert_matrix_type(mat_disk, type = "uint32_t")
  file_path = file.path(bpcells_dir, paste0("ref_", size_ref[size]))
  write_matrix_dir(
    mat = mat_disk,
    dir = file_path
  )
  mat <- open_matrix_dir(dir = file_path)
  data_ref <- CreateSeuratObject(counts = mat)
})

# Quality control ####
timers$with_timer("Quality control", {
  data_query[["percent.mt"]] <- PercentageFeatureSet(
    data_query, pattern = "^MT-")
  data_query <- subset(
    data_query, subset = nFeature_RNA > 200 & percent.mt < 5)
})

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
df = data_query@meta.data %>%
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
print(df, n = Inf)

write.csv(df, file.path(work_dir, "output", paste0(
  "test_transfer_seurat_bpcells_", size, "_accuracy.csv")), row.names = FALSE)

timers$print_summary(sort = FALSE)

timers_df <- timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$test <- 'test_transfer_seurat_bpcells'
timers_df$size <- size
write.csv(timers_df, output, row.names = FALSE)

unlink(file.path(bpcells_dir, size), recursive = TRUE)
unlink(file.path(bpcells_dir, paste0("ref_", size_ref[size])), recursive = TRUE)