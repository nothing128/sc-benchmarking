suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

work_dir <- 'output'

results <- bind_rows( # Use bind_rows() here
  lapply(
    list.files(path = "output", pattern = "test_", full.names = TRUE),
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

p <- ggplot(results, aes(x = duration, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_x_continuous(expand = c(0, 0)) +
  # coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Duration (s)", y = "", fill = "Number of cells") +
  theme_classic() +
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~size, scales = "free_x", ncol = 1)

ggsave("figures/comparison.png", plot = p, width = 10, height = 8)

my_colors <- c("test_basic_sc_multi_thread_subset" = "#ff0000",   
               "test_basic_sc_multi_thread_no_subset" = '#8b0000',
               "test_basic_scanpy_single_thread_subset" = "#CCCC00", 
               "test_seurat_bpcells_single_thread_subset" = "#0000ff",
               "test_basic_sc_single_thread_no_subset" = '#006400',
               "test_basic_sc_single_thread_subset" = '#00bb00',
               "test_basic_seurat_single_thread_subset" = "#00008B"
               
               ) 

p_1 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = duration, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(50, 100), scales = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(
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
  theme(legend.position = "right") +
  ggtitle("20K")

p_2 <- results %>% 
  filter(size == "400K") %>%
  ggplot(aes(x = duration, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(250, 400), scales = 0.5) +
  scale_x_break(c(1000, 2200), scales = 0.5) +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("400K")

p_3 <- results %>% 
  filter(size == "1M") %>%
  ggplot(aes(x = duration, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(450, 900), scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "fill_group") +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/comparison.png", plot = combined_plot, width = 10, height = 12)

p <- ggplot(results, aes(x = memory, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_x_continuous(expand = c(0, 0)) +
  # coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Duration (s)", y = "", fill = "Number of cells") +
  theme_classic() +
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~size, scales = "free_x", ncol = 1)

ggsave("figures/memory_comparison.png", plot = p, width = 10, height = 8)



my_colors <- c("test_basic_sc_multi_thread_subset" = "#ff0000",   
               "test_basic_sc_multi_thread_no_subset" = '#8b0000',
               "test_basic_scanpy_single_thread_subset" = "#CCCC00", 
               "test_seurat_bpcells_single_thread_subset" = "#0000ff",
               "test_basic_sc_single_thread_no_subset" = '#006400',
               "test_basic_sc_single_thread_subset" = '#00bb00',
               "test_basic_seurat_single_thread_subset" = "#00008B"
               
) 

p_1 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = memory, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  scale_fill_manual(
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
  theme(legend.position = "right") +
  ggtitle("20K")

p_2 <- results %>% 
  filter(size == "400K") %>%
  ggplot(aes(x = memory, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +

  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors) +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("400K")

p_3 <- results %>% 
  filter(size == "1M") %>%
  ggplot(aes(x = memory, y = operation, fill = fill_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "fill_group") +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/memory_comparison.png", plot = combined_plot, width = 10, height = 12)


