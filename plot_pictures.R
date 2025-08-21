suppressPackageStartupMessages({
  library(ggplot2)
  library(magick)
  library(ggimage)
  library(dplyr)
  library(purrr)
})

# --- Parameters ---
nn_val <- c(10,12,14)
enn_val <- c(26,28,30)
image_dir <- "figures"
crop_geometry <- "750x750+200+0" 

# --- Create a temporary directory for the cropped images ---
temp_crop_dir <- file.path(tempdir(), "cropped_images")
if (!dir.exists(temp_crop_dir)) {
  dir.create(temp_crop_dir)
}



plot_data <- expand.grid(enn = enn_val, nn = nn_val) %>%
  mutate(
    original_path = normalizePath(file.path(image_dir, paste0(nn, "neighbors_", enn, "extras.png"))),
    cropped_path = file.path(temp_crop_dir, basename(original_path))
  )

walk2(plot_data$original_path, plot_data$cropped_path, ~{
  image_read(.x) %>%
    image_crop(geometry = crop_geometry) %>%
    image_write(path = .y)
})



# --- Create the Plot ---

final_plot <- ggplot(plot_data, aes(x = 0, y = 0)) + 
  geom_image(aes(image = cropped_path), size = 1) +
  facet_grid(
    rows = vars(enn), 
    cols = vars(nn),
    labeller = labeller(
      nn = function(val) paste("nn:", val),
      enn = function(val) paste("Extras:", val)
    )
  ) +
  coord_fixed() +
  theme_void() +
  theme(
    panel.spacing = unit(0, "lines"), 
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 10, margin = margin(t = 3, b = 3)),
    plot.background = element_rect(color = "black", fill = NA)
  ) 

ggsave("final_output_3.png",final_plot)

