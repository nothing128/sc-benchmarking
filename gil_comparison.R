suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})
file_path <- 'output'

all_test_files <- list.files(
  path = file_path,
  pattern = "^test.*csv$", 
  full.names = TRUE
)

benchmarking_files <- grep(
  pattern = "_accuracy\\.csv$", 
  x = all_test_files,          
  invert = TRUE,               
  value = TRUE                 
)


benchmark_df <- bind_rows(
  lapply(
    benchmarking_files,
    function(file_path) {
      # Read CSV, keeping original names and strings as characters
      df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
      df$nogil <- if (grepl("nogil",file_path , fixed = TRUE))  "nogil" else "gil"
      
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
basic_levels <- c("Load data",
                  "Quality control",
                  "Doublet detection",
                  "Feature selection",
                  "Normalization",
                  "PCA",
                  "Nearest neighbors",
                  "Clustering (3 res.)",
                  "Embedding",
                  "Plot embedding",
                  "Find markers")
benchmark_df$subset <- ifelse(benchmark_df$subset == "true", "subset", "no subset")
benchmark_df <- benchmark_df %>% filter(subset=="subset")
p <- ggplot(benchmark_df, aes(x = duration, y = operation, fill = interaction(subset, num_threads, nogil))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_discrete(limits = rev(basic_levels)) +
  scale_x_continuous(expand = c(0, 0)) +
  # coord_cartesian(xlim = c(0, 50)) +
  labs(x = "Duration (s)", y = "", fill = "options") +
  theme_classic() +
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~size, scales = "free_x", ncol = 1)

ggsave("figures/comparison_gil_1.png", plot = p, width = 10, height = 8)
p
