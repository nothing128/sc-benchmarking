suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

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
results <- subset(results, select = -c(aborted,subset))



results <- results %>% mutate(
  thread_str = coalesce(as.character(num_threads), "single_thread"),
)
results$thread_str[results$thread_str == '1'] <- 'single_thread'
results$thread_str[results$thread_str == '-1'] <- 'multi_thread'
results$fill_group = paste(results$test, results$thread_str, sep = "_")
results$operation[results$operation == 'Load data (h5ad/rds)'] <- 'Load data'
results$operation[results$operation == 'Clustering (3 resolutions)'] <- 'Clustering (3 res.)'

results <- results %>%
  group_by(fill_group,size) %>%
  mutate(total = sum(duration)) %>%
  ungroup()
results<-subset(results, operation == "PCA")
results<- results[c("total","operation","size","fill_group")]
results$operation<-"total runtime"
results$duration<-results$total
results$total=NULL

results <- results %>% 
  mutate(
    size = 
      factor(
        size, 
        levels = c("20K",
                   "400K",
                   "1M"
                   )))

my_colors <- c("test_basic_sc_multi_thread" = "#aec6e7",   
               "test_basic_scanpy_single_thread" = "#bcbd23", 
               "test_seurat_bpcells_single_thread" = "#97df89",
               "test_basic_sc_single_thread" = '#ffbb77',
               "test_basic_seurat_single_thread" = "#2ba02d"
               
) 

plot<- results %>%
  ggplot(aes(x = size, y = duration, color=fill_group,group=fill_group))+
  geom_line() +
  theme_classic() +
  scale_color_manual(
    name = "model type",
    values = my_colors,
    labels = c("test_basic_sc_multi_thread" = "multithreaded brisc",   
               "test_basic_scanpy_single_thread" = "scanpy", 
               "test_seurat_bpcells_single_thread" = "seurat using bpcells",
               "test_basic_sc_single_thread" = 'single threaded brisc',
               "test_basic_seurat_single_thread" = "seurat"
               
    ) # Item labels
  ) +
  theme(legend.position = "right")
ggsave("figures/runtime_against_size.png", plot = plot, width = 10, height = 6)

