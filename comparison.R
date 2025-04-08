suppressPackageStartupMessages({
  library(tidyverse)
  library(ggbreak)
})

work_dir = 'projects/sc-benchmarking'

results = do.call(
    rbind, 
    lapply(list.files(path = paste0(work_dir, "/output"), 
    pattern = "test_basic_", full.names = TRUE), read.csv)) %>%
  mutate(test_size = paste(test, size, sep = "_"))

p <- ggplot(results, aes(x = duration, y = operation, fill = test_size)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_break(c(50000, 145000)) +
  labs(x = "Duration (ms)", y = "",
    fill = "Test & Number of cells") +
  theme_classic() +
  theme(
    axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(paste0(work_dir, "/figures/comparison.png"), 
    plot = p, width = 10, height = 8)