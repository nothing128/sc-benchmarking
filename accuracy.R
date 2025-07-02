suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
  library(aplot)
  library(cowplot)  
  library(dplyr)
})

data_dir <- "output4"

# 2. Get a list of all CSV file paths
files_to_read <- list.files(
  path = data_dir,
  pattern = "\\_accuracy.csv$",
  full.names = TRUE
)

# 3. Use lapply to read each file and add the new column
list_of_data_frames <- lapply(files_to_read, function(file_path) {
  # Read the CSV file
  df <- read.csv(file_path)
  
  # Add a new column with the filename (using basename to strip the directory path)
  df$source_file <- basename(file_path)
  
  # Return the modified data frame
  return(df)
})
my_colors <- c("test_transfer_sc_20K_multi_thread_subset_accuracy.csv"  = '#9bcae1',
               "test_transfer_sc_400K_multi_thread_subset_accuracy.csv" = "#6bafd6",
               "test_transfer_sc_1M_multi_thread_subset_accuracy.csv"   = "#3282bd",
               "test_transfer_scanpy_20K_accuracy.csv"                  = "#dadb8d",
               "test_transfer_scanpy_400K_accuracy.csv"                 = '#bcbd23',
               "test_transfer_scanpy_1M_accuracy.csv"                   = '#695a15',
               "test_transfer_seurat_20K_accuracy.csv"                  = '#a1d99a',
               "test_transfer_seurat_400K_accuracy.csv"                 = '#73c475',
               "test_transfer_seurat_1M_accuracy.csv"                   = '#30a355'
)
my_labels <- c(
  "test_transfer_sc_20K_multi_thread_subset_accuracy.csv" = "brisc",
  "test_transfer_sc_400K_multi_thread_subset_accuracy.csv" = "brisc",
  "test_transfer_sc_1M_multi_thread_subset_accuracy.csv" = "brisc",
  "test_transfer_scanpy_20K_accuracy.csv" = "Scanpy",
  "test_transfer_scanpy_400K_accuracy.csv" = "Scanpy",
  "test_transfer_scanpy_1M_accuracy.csv" = "Scanpy",
  "test_transfer_seurat_20K_accuracy.csv" = "Seurat",                
  "test_transfer_seurat_400K_accuracy.csv" = "Seurat",                
  "test_transfer_seurat_1M_accuracy.csv" = "Seurat"
)
# 4. Combine the list of data frames into a single data frame
data <- do.call(rbind, list_of_data_frames)
p_1 <- data %>% 
  filter(grepl("20K",source_file)) %>% 
  filter(!grepl("bpcells",source_file)) %>%
  ggplot(aes(x = subclass, y = percent_correct, title="20K", fill=source_file)) +
  ylab("correct %")+
  labs(title="20K dataset")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(
    name = "toolkit",
    labels= my_labels,
    values = my_colors)+
  theme_classic() +
  theme(legend.position = "bottom")
p_2 <- data %>% 
  filter(grepl("400K",source_file)) %>%
  filter(!grepl("bpcells",source_file)) %>%
  ggplot(aes(x = subclass, y = percent_correct, fill=source_file)) +
  ylab("correct %")+
  labs(title="400K dataset")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(
    name = "toolkit",
    labels= my_labels,
    values = my_colors)+
  theme_classic() +
  theme(legend.position = "bottom")
p_3 <- data %>% 
  filter(grepl("1M",source_file)) %>%
  filter(!grepl("bpcells",source_file)) %>%
  ggplot(aes(x = subclass, y = percent_correct, fill=source_file)) +
  ylab("correct %")+
  labs(title="1M dataset")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(
    name = "toolkit",
    labels= my_labels,
    values = my_colors)+
  theme_classic() +
  theme(legend.position = "bottom")
combined_plot <- aplot::plot_list(p_1, p_2, p_3, ncol = 1)
ggsave("figures/accuracy.png", plot = combined_plot, width = 20, height = 10)

