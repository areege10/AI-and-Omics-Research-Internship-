# -------------------------------
# Define classification function
# -------------------------------
classify_gene <- function(logFC, padj) {
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# Example test calls
classify_gene(logFC = 2, padj = 0.01)
classify_gene(logFC = -2, padj = 0.01)
classify_gene(logFC = 0.5, padj = 0.5)

# -------------------------------
# Input and output directories
# -------------------------------
input_dir <- "raw_data"
output_dir <- "Results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# -------------------------------
# Files to process
# -------------------------------
files_to_process <- c("DEGs_data_1.csv", "DEGs_data_2.csv")
result_list <- list()

for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  
  input_file_path <- file.path(input_dir, file_name)
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing padj values...\n")
  
  # Replace missing padj with 1
  missing_count <- sum(is.na(data$padj))
  cat("Missing values in 'padj':", missing_count, "\n")
  data$padj[is.na(data$padj)] <- 1
  
  # Apply classification
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  cat("Gene classification completed.\n")
  
  # Save results in R
  result_list[[file_name]] <- data
  
  # Save results in Results folder
  output_file_path <- file.path(output_dir, paste0("DEGs_results_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  # Print summary
  cat("Summary counts:\n")
  print(table(data$status))
}

# -------------------------------
# Access results back from list
# -------------------------------
results_1 <- result_list[[1]]
results_2 <- result_list[[2]]


save.image(file = "YourName_Class_2_Assignment.RData")
