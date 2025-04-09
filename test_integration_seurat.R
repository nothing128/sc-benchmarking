suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})
options(Seurat.object.assay.version = "v5")

work_dir = "projects/sc-benchmarking"
data_dir = "scratch/single-cell/SEAAD"
source(file.path(work_dir, "utils_local.R"))

system_info()

size = "20K"

timers = TimerCollection(silent = FALSE)

timers$with_timer("Load data", {
    data <- readRDS(paste0(data_dir, "/SEAAD_raw_", size, ".rds"))
})

timers$with_timer('Split data', {
    data[["RNA"]] <- split(data[['RNA']], f = data$ad_dx)
})
