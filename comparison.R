suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
})

work_dir <- "projects/sc-benchmarking"

results <- bind_rows(
  lapply(
    list.files(path = paste0(work_dir, "/output"),
    pattern = "\\.csv$", full.names = TRUE),
    function(file_path) {
      read.csv(file_path, stringsAsFactors = FALSE)
    }
  )
)

results <- results %>% mutate(
  thread_str = ifelse(is.na(num_threads), "single_thread",
    ifelse(num_threads == 1, "single_thread", "multi_thread")),
  subset_str = ifelse(is.na(subset), "subset",
    ifelse(subset == TRUE, "subset", "no_subset")),
) %>%
  mutate(
    fill_group = paste(test, thread_str, subset_str, sep = "_")
  )

results <- results %>% 
  mutate(operation = factor(operation, 
    levels = c("Load data",
              "Load data query",
              "Load data ref",
              "Quality control",
              "Doublet detection",
              "Feature selection",
              "Normalization",
              "PCA",
              "Neighbor graph",
              "Clustering (3 res.)",
              "Embedding",
              "Plot embeddings",
              "Find markers",
              "Align datasets",
              "Transfer labels"
  )))

basic_results <- results %>% filter(grepl("basic", test))
transfer_results <- results %>% filter(grepl("transfer", test))

my_colors <- c(
  "test_basic_sc_multi_thread_subset" = "#aec6e7",   
  "test_basic_sc_multi_thread_no_subset" = '#1f78b4',
  "test_basic_scanpy_single_thread_subset" = "#bcbd23", 
  "test_seurat_bpcells_single_thread_subset" = "#97df89",
  "test_basic_sc_single_thread_no_subset" = '#fe7f0e',
  "test_basic_sc_single_thread_subset" = '#ffbb77',
  "test_basic_seurat_single_thread_subset" = "#2ba02d",
  "test_transfer_sc_multi_thread_subset" = "#d62728",
  "test_transfer_sc_multi_thread_no_subset" = "#9467bd",
  "test_transfer_sc_single_thread_subset" = "#8c564b",
  "test_transfer_sc_single_thread_no_subset" = "#e377c2",
  "test_transfer_scanpy_single_thread_subset" = "#7f7f7f",
  "test_transfer_seurat_single_thread_subset" = "#bcbd22"
  ) 

# runtime plots
# split by dataset size
p_1 <- basic_results %>% 
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

p_2 <- basic_results %>% 
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

p_3 <- basic_results %>% 
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
ggsave("figures/comparison_basic.png", plot = combined_plot, width = 10, height = 12)

# memory comparison plot
p_1 <- basic_results %>% 
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

p_2 <- basic_results %>% 
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

p_3 <- basic_results %>% 
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
ggsave("figures/memory_comparison_basic.png", plot = combined_plot, width = 10, height = 12)


# Transfer plots
p_1 <- transfer_results %>%
    filter(size == "20K") %>%
    ggplot(aes(x = duration, y = operation, fill = fill_group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(results$operation))) +
    scale_fill_manual(
        name = "model type",
        values = my_colors,
        labels = c("test_transfer_sc_multi_thread_subset" = "sc multi-thread subset",
                   "test_transfer_sc_multi_thread_no_subset" = "sc multi-thread no-subset",
                   "test_transfer_sc_single_thread_subset" = "sc single-thread subset",
                   "test_transfer_sc_single_thread_no_subset" = "sc single-thread no-subset",
                   "test_transfer_scanpy_single_thread_subset" = "scanpy",
                   "test_transfer_seurat_single_thread_subset" = "seurat"
                   )
    ) +
    theme_classic() +
    theme(legend.position = "right") +
    ggtitle("20K")

p_2 <- transfer_results %>%
    filter(size == "400K") %>%
    ggplot(aes(x = duration, y = operation, fill = fill_group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(results$operation))) +
    scale_fill_manual(values = my_colors) +
    labs(x = "Duration (s)", y = "") +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("400K")

p_3 <- transfer_results %>%
    filter(size == "1M") %>%
    ggplot(aes(x = duration, y = operation, fill = fill_group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(results$operation))) +
    scale_fill_manual(values = my_colors, name = "fill_group") +
    labs(x = "Duration (s)", y = "") +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/comparison_transfer.png", plot = combined_plot, width = 10, height = 12)


p_1 <- transfer_results %>%
    filter(size == "20K") %>%
    ggplot(aes(x = memory, y = operation, fill = fill_group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(levels(results$operation))) +
    labs(title="memory usage", x="memory usage (GiB)", y="") +
    scale_fill_manual(
        name = "model type",
        values = my_colors,
        labels = c("test_transfer_sc_multi_thread_subset" = "sc multi-thread subset",
                   "test_transfer_sc_multi_thread_no_subset" = "sc multi-thread no-subset",
                   "test_transfer_sc_single_thread_subset" = "sc single-thread subset",
                   "test_transfer_sc_single_thread_no_subset" = "sc single-thread no-subset",
                   "test_transfer_scanpy_single_thread_subset" = "scanpy",
                   "test_transfer_seurat_single_thread_subset" = "seurat"
                   )
    ) +
    theme_classic() +
    theme(legend.position = "right") +
    ggtitle("20K")

p_2 <- transfer_results %>%
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

p_3 <- transfer_results %>%
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
ggsave("figures/memory_comparison_transfer.png", plot = combined_plot, width = 10, height = 12)


