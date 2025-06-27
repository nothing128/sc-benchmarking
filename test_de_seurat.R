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
timers = TimerMemoryCollection(silent = TRUE)

# Note: Loading is much slower from $SCRATCH disk

# Load data ####
timers$with_timer("Load data", {
  data <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
})

# Not timed
data[["nCount_RNA"]] <- colSums(data@assays$RNA@counts)
data[["nFeature_RNA"]] <- colSums(data@assays$RNA@counts > 0)

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
timers_df$test <- 'test_basic_seurat'
timers_df$size <- size

write.csv(timers_df, output, row.names = FALSE)
'''
--- Timing Summary ---
Load data took 11s 755ms (1.3%) using 0.00 GiB (1.4%)
Quality control took 5s 785ms (0.6%) using 0.00 GiB (2.6%)
Normalization took 7s 818ms (0.8%) using 0.00 GiB (2.9%)
Differential expression (wilcoxon) took 34s 453ms (3.7%) using 0.00 GiB (3.1%)
Pseudobulk took 8s 14ms (0.9%) using 0.00 GiB (3.5%)
Differential expression (DESeq2) took 14m 20s (92.7%) using 0.00 GiB (2.8%)

Total time: 15m 28s
'''