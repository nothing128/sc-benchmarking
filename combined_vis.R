suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

file_path <- 'poster_output'

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
my_colors2 <- c("brisc single-threaded"= "#ffbb77",
                "brisc multi-threaded" = "#aec6e7",
                "scanpy"="#bcbd23",
                "seurat"= "#97df89")

benchmark_df$operation[benchmark_df$operation == 'Data transformation (pseudobulk / normalization)'] <- 'Data transformation'
benchmark_df <- benchmark_df %>%
  unite(
    group,            
    library, num_threads,  
    sep = " ",
    remove = FALSE,       
    na.rm = TRUE
  )


de_levels <- c("Load data",
               "Quality control",
               "Doublet detection",
               "Differential expression",
               "Data transformation")

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

transfer_levels<- c("Load data (query)",
                     "Load data (ref)",
                     "Quality control",
                     "Normalization",
                     "Feature selection",
                     "PCA",
                     "Transfer labels")

levels <- list("de"=de_levels,
            "basic"=basic_levels,
            "transfer"=transfer_levels)

create_plot <- function(df, size_arg, test_arg, runtime, title=size_arg, breaks=NULL,scales=NULL) {
  break_layer <- if (!is.null(breaks)) {
    scale_x_break(breaks = breaks, scales = scales)
  } else {
    NULL
  }
  legend_theme <- if (size_arg == "20K") {
    theme(legend.position = "right")
  } else {
    theme(legend.position = "none")
  }
  plot <- df %>%
    filter(size == .env$size_arg) %>%
    filter(test == .env$test_arg) %>%
    mutate(
      operation = 
        factor(
          operation, 
          levels = levels[[test_arg]])) %>%
    ggplot(aes(x = .data[[if (runtime) "duration" else "memory"]],
               y = operation, fill=group)) + 
    geom_bar(stat = "identity", position = "dodge") +
    break_layer +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = my_colors2) +
    scale_y_discrete(limits = rev(levels[[test_arg]])) +
    theme_classic() +
    legend_theme+
    labs(fill="package")
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}
p_2 <- create_plot(benchmark_df,"400K","basic",TRUE, "brisc is significantly faster than scanpy and seurat")
p_2
p_2 <- create_plot(benchmark_df,"400K","basic",TRUE, "", c(150, 400,650,1000),c(0.5,0.5))
p_2
p_2 <- create_plot(benchmark_df,"400K","basic",FALSE, "")
p_2




p_1 <- create_plot(benchmark_df,"20K","de",TRUE,"de 20K")
p_2 <- create_plot(benchmark_df,"20K","transfer",TRUE,"transfer 20K",c(120,999),c(0))
p_3 <- create_plot(benchmark_df,"20K","basic",TRUE,"basic 20K",c(100,200),c("free"))
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/20K_runtime.png",plot = combined_plot, width = 10, height = 12)


# I don't like this but each plot has a unique break
p_1 <- create_plot(benchmark_df,"20K","de",TRUE,c(50, 600),c(0.5))
p_2 <- create_plot(benchmark_df,"400K","de",TRUE,c(50, 700),c(0.5))
p_3 <- create_plot(benchmark_df,"1M","de",TRUE,c(500, 3000),c("free"))
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/de_runtime.png",plot = combined_plot, width = 10, height = 12)

p_1 <- create_plot(benchmark_df,"20K","de",FALSE)
p_2 <- create_plot(benchmark_df,"400K","de",FALSE)
p_3 <- create_plot(benchmark_df,"1M","de",FALSE)
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/de_memory.png",plot = combined_plot, width = 10, height = 12)

p_1 <- create_plot(benchmark_df,"20K","transfer",TRUE)
p_2 <- create_plot(benchmark_df,"400K","transfer",TRUE,c(50, 900),c(0.5))
p_3 <- create_plot(benchmark_df,"1M","transfer",TRUE,c(220, 500),c(0.5))
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/transfer_runtime.png",plot = combined_plot, width = 10, height = 12)

p_1 <- create_plot(benchmark_df,"20K","transfer",FALSE)
p_2 <- create_plot(benchmark_df,"400K","transfer",FALSE)
p_3 <- create_plot(benchmark_df,"1M","transfer",FALSE)
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/transfer_memory.png",plot = combined_plot, width = 10, height = 12)

p_1 <- create_plot(benchmark_df,"20K","basic",TRUE,c(50, 100),c(0.5))
p_2 <- create_plot(benchmark_df,"400K","basic",TRUE,c(250, 400,1000,2200),c(0.5,0.5))
p_3 <- create_plot(benchmark_df,"1M","basic",TRUE,c(450,900),c("free"))
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/basic_runtime.png",plot = combined_plot, width = 10, height = 12)

p_1 <- create_plot(benchmark_df,"20K","basic",FALSE)
p_2 <- create_plot(benchmark_df,"400K","basic",FALSE)
p_3 <- create_plot(benchmark_df,"1M","basic",FALSE)
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/basic_memory.png",plot = combined_plot, width = 10, height = 12)
  

