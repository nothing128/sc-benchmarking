# test_read_obs.R
# This script is for testing the read_h5ad_obs function in utils_local.R

# Load utilities
source("projects/sc-benchmarking/utils_local.R")

# Define path to test data
# Using a small dataset for quick testing.
# This assumes the script is run from the workspace root.
data_dir <- "single-cell/SEAAD"
size <- "20K"
h5ad_path <- file.path(data_dir, paste0("SEAAD_raw_", size, ".h5ad"))

cat("--- Testing read_h5ad_obs function ---\n")
cat("Using file:", h5ad_path, "\n\n")

# Check if file exists
if (!file.exists(h5ad_path)) {
  stop("Test file not found at:", h5ad_path)
}

# Run the function and capture potential errors
obs_metadata <- NULL
tryCatch({
  obs_metadata <- read_h5ad_obs(h5ad_path)
}, error = function(e) {
  cat("An error occurred while running read_h5ad_obs:\n")
  print(e)
})

# If the function ran successfully, print diagnostics
if (!is.null(obs_metadata)) {
  cat("\n--- Function executed successfully ---\n")
  cat("\nStructure of the returned data.frame (obs_metadata):\n")
  str(obs_metadata)
  
  cat("\n\nFirst 6 rows of the data.frame:\n")
  print(head(obs_metadata))
  
  cat("\n\nColumn names:\n")
  print(colnames(obs_metadata))
  
  # Check for any column names still containing "/codes"
  if (any(grepl("/codes", colnames(obs_metadata)))) {
    cat("\n\n[WARNING] Found column names containing '/codes'. The function did not clean them up correctly.\n")
  } else {
    cat("\n\n[SUCCESS] No column names with '/codes' found.\n")
  }
} else {
  cat("\n--- Function execution failed ---\n")
}

cat("\n--- Test script finished ---\n") 