suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})
options(Seurat.object.assay.version = "v5")

work_dir = "projects/rrg-wainberg/karbabi/sc-benchmarking"
data_dir = "single-cell/SEAAD"
source(file.path(work_dir, "utils_local.R"))

system_info()

size = "400K"

timers = TimerCollection(silent = FALSE)

timers$with_timer("Load query data", {
    data <- readRDS(
        paste0(data_dir, "/subsampled/SEAAD_raw_", size, ".rds"))
})

timers$with_timer("Load reference data", {
    ref <- readRDS(
        paste0(data_dir, "/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"))
})

timers$with_timer('Split data', {
    data[["RNA"]] <- split(data[['RNA']], f = data$ad_dx)
})
