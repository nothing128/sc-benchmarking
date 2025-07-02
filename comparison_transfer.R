suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

work_dir <- 'output4'

all_test_files <- list.files(
  path = work_dir,
  pattern = "^test_transfer.*\\.csv$", # A simple pattern list.files can handle
  full.names = TRUE
)

final_file_list <- grep(
  pattern = "_accuracy\\.csv$", # The pattern we want to exclude
  x = all_test_files,          # The list of files to search through
  invert = TRUE,               # Invert the match: return non-matches
  value = TRUE                 # Return the file names themselves
)


# --- The rest of your code remains the same ---
# It now uses the correctly filtered 'final_file_list'

results <- bind_rows(
  lapply(
    final_file_list, # Use the correctly filtered list
    function(file_path) {
      # Read CSV, keeping original names and strings as characters
      df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
      
      # --- Standardize column names (stripping quotes) ---
      current_names <- names(df)
      new_names <- sapply(current_names, function(col_name) {
        if (startsWith(col_name, "\"") && endsWith(col_name, "\"") && nchar(col_name) > 1) {
          return(substr(col_name, 2, nchar(col_name) - 1))
        } else {
          return(col_name)
        }
      })
      names(df) <- new_names
      
      aborted_col_standardized_name <- "aborted"
      
      if (aborted_col_standardized_name %in% names(df)) {
        df[[aborted_col_standardized_name]] <- as.character(df[[aborted_col_standardized_name]])
      }
      return(df)
    }
  )
)
results <- subset(results, select = -c(aborted))
results <- subset(results, select = -c(num_threads))
results <- subset(results, select = -c(subset))
results <- results %>%
  filter(!str_detect(operation, "Doublet detection"))
results <- results %>% 
  mutate(
    operation = 
      factor(
        operation, 
        levels = c("Load data (query)",
                   "Load data (ref)",
                   "Quality control",
                   "Normalization",
                   "Feature selection",
                   "PCA",
                   "Transfer labels")))




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

ggsave("figures/comparison_transfer_time.png", plot = p, width = 10, height = 8)

my_colors<- c("test_transfer_sc" = "#aec6e7",   
              "test_transfer_scanpy" = "#bcbd23", 
              "test_transfer_seurat" = "#2ba02d",
              "test_transfer_seurat_bpcells" = "#97df89")

p_1 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  labs(x = "Duration (s)", y = "") +
  scale_fill_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_transfer_sc" = "brisc",   
               "test_transfer_scanpy" = "scanpy", 
               "test_transfer_seurat" = "seurat",
               "test_transfer_seurat_bpcells" = "seurat with bpcells") # Item labels
  ) +
  theme_classic() +
  theme(legend.position = "right") +
  ggtitle("20K")

p_2 <- results %>% 
  filter(size == "400K") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_break(c(500, 900), scales = 0.5) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("400K")

p_3 <- results %>% 
  filter(size == "1M") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_break(c(220, 500), scales = 0.5) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "test") +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/comparison_transfer_time.png", plot = combined_plot, width = 10, height = 12)

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

ggsave("figures/memory_comparison_transfer.png", plot = p, width = 10, height = 8)


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
    labels = c("test_transfer_sc" = "brisc",   
               "test_transfer_scanpy" = "scanpy", 
               "test_transfer_seurat" = "seurat",
               "test_transfer_seurat_bpcells" = "seurat with bpcells") # Item labels
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
  scale_fill_manual(values = my_colors) +
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
ggsave("figures/memory_comparison_transfer.png", plot = combined_plot, width = 10, height = 12)


