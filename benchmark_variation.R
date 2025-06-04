suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

work_dir <- 'output2'

results <- bind_rows( # Use bind_rows() here
  lapply(
    list.files(path = work_dir, pattern = "test_", full.names = TRUE),
    function(file_path) {
      # Read CSV, keeping original names and strings as characters
      df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
      
      # --- Step 1: Standardize column names (stripping quotes) ---
      current_names <- names(df)
      new_names <- sapply(current_names, function(col_name) {
        if (startsWith(col_name, "\"") && endsWith(col_name, "\"") && nchar(col_name) > 1) {
          return(substr(col_name, 2, nchar(col_name) - 1))
        } else {
          return(col_name)
        }
      })
      names(df) <- new_names # Standardized names are now active in df
      aborted_col_standardized_name <- "aborted"
      
      if (aborted_col_standardized_name %in% names(df)) {
        # Convert the 'aborted' column to character in this specific data frame
        df[[aborted_col_standardized_name]] <- as.character(df[[aborted_col_standardized_name]])
      }
      return(df)
    }
  )
) 
results <- subset(results, select = -c(aborted))



results <- results %>% mutate(
  thread_str = coalesce(as.character(num_threads), "single_thread"),
  subset_str = coalesce(as.character(subset), "subset"),
)
results$thread_str[results$thread_str == '1'] <- 'single_thread'
results$thread_str[results$thread_str == '-1'] <- 'multi_thread'
results$subset_str[results$subset_str == 'false'] <- 'no_subset'
results$subset_str[results$subset_str == 'true'] <- 'subset'
results$fill_group = paste(results$test, results$thread_str, results$subset_str, sep = "_")
results$operation[results$operation == 'Load data (h5ad/rds)'] <- 'Load data'
results$operation[results$operation == 'Clustering (3 resolutions)'] <- 'Clustering (3 res.)'

results <- results %>% 
  mutate(
    operation = 
      factor(
        operation, 
        levels = c("Load data",
                   "Quality control",
                   "Doublet detection",
                   "Feature selection",
                   "Normalization",
                   "PCA",
                   "Neighbor graph",
                   "Clustering (3 res.)",
                   "Embedding",
                   "Plot embeddings",
                   "Find markers")))

my_colors <- c("test_basic_sc_multi_thread_subset" = "#aec6e7",   
               "test_basic_sc_multi_thread_no_subset" = '#1f78b4',
               "test_basic_scanpy_single_thread_subset" = "#bcbd23", 
               "test_seurat_bpcells_single_thread_subset" = "#97df89",
               "test_basic_sc_single_thread_no_subset" = '#fe7f0e',
               "test_basic_sc_single_thread_subset" = '#ffbb77',
               "test_basic_seurat_single_thread_subset" = "#2ba02d"
               
)

p_1 <- results %>% 
  filter((test != "test_basic_scanpy") & (test != "test_seurat_bpcells")) %>%
  ggplot(aes(x = duration, y = operation, color = fill_group))+
  geom_point(stat = "identity") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_break(c(75, 200), scales = 0.5) +
  coord_cartesian(xlim = c(0, 250)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_color_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc_multi_thread_subset" = "multithreaded with subset",   
               "test_basic_sc_multi_thread_no_subset" = "multithreaded without subset",
               "test_basic_scanpy_single_thread_subset" = "scanpy", 
               "test_seurat_bpcells_single_thread_subset" = "seurat using bpcells",
               "test_basic_sc_single_thread_no_subset" = 'single threaded without subset',
               "test_basic_sc_single_thread_subset" = 'single threaded with subset',
               "test_basic_seurat_single_thread_subset" = "seurat"
               
    ) # Item labels
  ) +
  theme_classic() +
  theme(legend.position = "right")

p_2 <- results %>% 
  filter((test == "test_basic_scanpy") | (test == "test_seurat_bpcells")) %>%
  ggplot(aes(x = duration, y = operation, color = fill_group))+
  geom_point(stat = "identity") + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_break(c(200, 400), scales = 0.5) +
  scale_x_break(c(1200, 2300), scales = 0.5) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_color_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc_multi_thread_subset" = "multithreaded with subset",   
               "test_basic_sc_multi_thread_no_subset" = "multithreaded without subset",
               "test_basic_scanpy_single_thread_subset" = "scanpy", 
               "test_seurat_bpcells_single_thread_subset" = "seurat using bpcells",
               "test_basic_sc_single_thread_no_subset" = 'single threaded without subset',
               "test_basic_sc_single_thread_subset" = 'single threaded with subset',
               "test_basic_seurat_single_thread_subset" = "seurat"
               
    ) # Item labels
  ) +
  theme_classic() +
  theme(legend.position = "right")
combined_plot <- aplot::plot_list(p_1, p_2, ncol = 1)
ggsave("figures/runtime_variation.png", plot = combined_plot, width = 10, height = 6)


p_1_mem <- results %>% 
  filter((test != "test_basic_scanpy") & (test != "test_seurat_bpcells")) %>%
  ggplot(aes(x = memory, y = operation, color = fill_group))+
  geom_point(stat = "identity") + 
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 50)) +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_color_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc_multi_thread_subset" = "multithreaded with subset",   
               "test_basic_sc_multi_thread_no_subset" = "multithreaded without subset",
               "test_basic_scanpy_single_thread_subset" = "scanpy", 
               "test_seurat_bpcells_single_thread_subset" = "seurat using bpcells",
               "test_basic_sc_single_thread_no_subset" = 'single threaded without subset',
               "test_basic_sc_single_thread_subset" = 'single threaded with subset',
               "test_basic_seurat_single_thread_subset" = "seurat"
               
    ) # Item labels
  ) +
  theme_classic() +
  theme(legend.position = "right")

p_2_mem <- results %>% 
  filter((test == "test_basic_scanpy") | (test == "test_seurat_bpcells")) %>%
  ggplot(aes(x = memory, y = operation, color = fill_group))+
  geom_point(stat = "identity") + 
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_color_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc_multi_thread_subset" = "multithreaded with subset",   
               "test_basic_sc_multi_thread_no_subset" = "multithreaded without subset",
               "test_basic_scanpy_single_thread_subset" = "scanpy", 
               "test_seurat_bpcells_single_thread_subset" = "seurat using bpcells",
               "test_basic_sc_single_thread_no_subset" = 'single threaded without subset',
               "test_basic_sc_single_thread_subset" = 'single threaded with subset',
               "test_basic_seurat_single_thread_subset" = "seurat"
               
    ) # Item labels
  ) +
  theme_classic() +
  theme(legend.position = "right")
combined_plot <- aplot::plot_list(p_1_mem, p_2_mem, ncol = 1)
ggsave("figures/memory_variation.png", plot = combined_plot, width = 10, height = 6)
