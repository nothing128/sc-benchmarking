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

    file_path = paste0(bpcells_dir, size)
    write_matrix_dir(
      mat = mat_disk,
      dir = file_path
    )
    mat <- open_matrix_dir(dir = file_path)
    data <- CreateSeuratObject(counts = mat)
  })

data_tmp <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
data <- AddMetaData(
  object = data, metadata = data_tmp@meta.data[,
    !colnames(data_tmp@meta.data) %in% colnames(data@meta.data)])
rm(data_tmp); gc()

# Not timed
data <- AddMetaData(data, metadata = data.frame(
  nCount_RNA = colSums(data@assays$RNA@layers$counts),
  nFeature_RNA = colSums(data@assays$RNA@layers$counts > 0)
))

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

# Not timed
data$ad_dx <- ifelse(data$ad_dx == "1", "AD", "Control")

# Differential expression ####
timers$with_timer("Differential expression (wilcoxon)", {
  data$group <- paste(data$subclass, data$ad_dx, sep = "_")
  Idents(data) <- "group"

  de_list <- list()
  for (subclass in unique(data$subclass)) {
    de_list[[subclass]] <- FindMarkers(
      data, 
      ident.1 = paste(subclass, "AD", sep = "_"), 
      ident.2 = paste(subclass, "Control", sep = "_")
    )
  }
  de <- do.call(rbind, de_list)
})

# Pseudobulk ####
timers$with_timer("Pseudobulk", {
  data <- AggregateExpression(
    data, 
    assays = "RNA", 
    return.seurat = TRUE, 
    group.by = c("ad_dx", "sample", "subclass")
  )
})

# Differential expression ####
timers$with_timer("Differential expression (DESeq2)", {
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
timers_df <- timers$to_dataframe(unit = "s", sort = FALSE)
timers_df$test <- 'test_de_seurat_bpcells'
timers_df$size <- size

write.csv(timers_df, output, row.names = FALSE)
