suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(patchwork)
})

work_dir = "projects/def-wainberg/karbabi/sc-benchmarking"
source(file.path(work_dir, "utils_local.R"))

system_info()

timers = TimerCollection(silent = FALSE)

data <- Read10X(
    data.dir = "single-cell/SEAAD/subsampled/SEAAD_raw_20K",
    gene.column=1)
data <- CreateSeuratObject(counts = data)

data = Read10X_h5(
    filename = "single-cell/SEAAD/subsampled/SEAAD_raw_20k.h5")
data <- CreateSeuratObject(counts = data)








# create a timer collection
timers = TimerCollection(silent = FALSE)

# time some operations
timers$with_timer('load data', {
  df = data.frame(x = 1:1000000, y = rnorm(1000000))
})

timers$with_timer('calculate mean', {
  mean_val = mean(df$y)
})

timers$with_timer('filter values', {
  filtered = df[df$y > 0, ]
})

# print timing summary
timers$print_summary()

# get results as a dataframe
results_df = timers$to_dataframe()
print(results_df)