suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

work_dir <- 'output4'

results <- bind_rows( # Use bind_rows() here
  lapply(
    list.files(path = work_dir, pattern = "^test_de.*[a-zA-Z].csv$", full.names = TRUE),
    function(file_path) {
      print(file_path)
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
results$operation[results$operation=="Differential expression (wilcoxon)"] <- "Differential expression"
results$operation[results$operation=="Differential expression (DESeq2)"] <- "Differential expression"
results$operation[results$operation=="Quality control"] <- "Quality control (single cell)"
results <- results %>% 
  mutate(
    operation = 
      factor(
        operation, 
        levels = c("Load data",
                   "Quality control (single cell)",
                   "Doublet detection",
                   "Normalization",
                   "Pseudobulk",
                   "Quality control (pseudobulk)",
                   "Differential expression")))

p <- ggplot(results, aes(x = duration, y = operation, fill = test)) +
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

ggsave("figures/comparison_de.png", plot = p, width = 10, height = 8)

my_colors <- c("test_basic_sc" = "#aec6e7",   
               "test_basic_scanpy" = "#bcbd23", 
               "test_basic_seurat_bpcells_single_thread_subset" = "#97df89",
               "test_basic_seurat" = "#2ba02d")
               
 

p_1 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(50, 600), scales = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc" = "brisc",   
               "test_basic_scanpy"= "scanpy", 
               "test_basic_seurat_bpcells" = "seurat using bpcells",
               "test_basic_sc_single_thread_subset" = 'single threaded brisc',
               "test_basic_seurat" = "seurat"
               
    ) # Item labels
  ) +
  theme_classic() +
  theme(legend.position = "right") +
  ggtitle("20K")

p_2 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(50, 700), scales = 0.5) +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "fill_group") +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("400K")

p_3 <- results %>% 
  filter(size == "1M") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(500, 3000), scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "fill_group") +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/comparison_de.png", plot = combined_plot, width = 10, height = 12)

results$memory <- ifelse(
  results$memory_unit == "KiB",               
  results$memory / (1024 * 1024),             
  results$memory                              
)

p <- ggplot(results, aes(x = memory, y = operation, fill = test)) +
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

ggsave("figures/memory_comparison_de.png", plot = p, width = 10, height = 8)


p_1 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = memory, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  scale_fill_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc" = "brisc",   
               "test_basic_scanpy"= "scanpy", 
               "test_basic_seurat_bpcells" = "seurat using bpcells",
               "test_basic_sc_single_thread_subset" = 'single threaded brisc',
               "test_basic_seurat" = "seurat"
               
    )
  ) +
  theme_classic() +
  theme(legend.position = "right") +
  ggtitle("20K")

p_2 <- results %>% 
  filter(size == "400K") %>%
  ggplot(aes(x = memory, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "test") +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("400K")

p_3 <- results %>% 
  filter(size == "1M") %>%
  ggplot(aes(x = memory, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "test") +
  labs(title="memory usage", x="memory usage (GiB)", y="") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/memory_comparison_de.png", plot = combined_plot, width = 10, height = 12)


