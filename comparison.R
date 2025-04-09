suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
})

work_dir <- 'scratch/projects/sc-benchmarking'
setwd(work_dir)

results <- do.call(
    rbind, 
    lapply(list.files(path = "output", 
    pattern = "test_basic_", full.names = TRUE), read.csv)) %>%
    filter(!operation %in% c("Load data (10X mtx)", "Load data (h5)")) %>%
    mutate(size = factor(size, levels = c("20K", "400K", "1M")),
           operation = factor(operation, 
            levels = c("Load data (h5ad/rds)",
              "Quality control",
              "Doublet detection",
              "Feature selection",
              "Normalization",
              "PCA",
              "Neighbor graph",
              "Clustering (3 resolutions)",
              "Embedding",
              "Plot embeddings",
              "Find markers")))

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

ggsave("figures/comparison.png", plot = p, width = 10, height = 8)

my_colors <- c("test_basic_sc" = "#FF6B6B",    
               "test_basic_scanpy" = "#4ECDC4", 
               "test_basic_seurat" = "#45B3E0",
               'test_basic_sc_single_thread' = '#FFD700') 

p_1 <- results %>% 
  filter(size == "20K") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(50, 100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "right") +
  ggtitle("20K")

p_2 <- results %>% 
  filter(size == "400K") %>%
  ggplot(aes(x = duration, y = operation, fill = test)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(200, 2500), scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
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
  scale_x_break(c(30, 40), scales = "free") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev(levels(results$operation))) +
  scale_fill_manual(values = my_colors, name = "test") +
  labs(x = "Duration (s)", y = "") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle("1M")

combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/comparison.png", plot = combined_plot, width = 10, height = 12)


