# SCRIPT TO CREATE EXCEL FILES FROM plot_data.csv
# Save this code in a separate file, e.g., create_excel_tables.R

# Load required libraries
if (!require("writexl")) {
  install.packages("writexl")
  library(writexl)
}

if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require("tidyr")) {
  install.packages("tidyr")
  library(tidyr)
}

# Read data from plot_data.csv
plot_data <- read.csv("D:/AcMet/3AcMet/epigenetic-analysis/Methylation/results/plot_data.csv", stringsAsFactors = FALSE)

cat("[INFO] Data read from plot_data.csv\n")
cat("Dataset dimensions:", dim(plot_data), "\n")
cat("Columns:", names(plot_data), "\n")

# Check unique values for Gene and Project
cat("Genes present:", unique(plot_data$Gene), "\n")
cat("Projects present:", unique(plot_data$Project), "\n")

# Define the desired order for the project columns (primary tumors)
project_col_order <- c("TCGA-LUAD", "TCGA-PRAD", "TCGA-COAD")
cat("\n[INFO] Desired project column order:", project_col_order, "\n")

# Create separate tables for CGAS and STING in wide format
create_gene_tables <- function(plot_data, order_vec) {
  # Filter and transform data for CGAS
  cgas_table <- plot_data %>%
    filter(Gene == "CGAS") %>%
    select(Probe, Project, BetaValue) %>%
    pivot_wider(
      names_from = Project,
      values_from = BetaValue,
      values_fn = mean  # Use mean if there are duplicates
    ) %>%
    select(Probe, all_of(order_vec))
  
  # Filter and transform data for STING
  sting_table <- plot_data %>%
    filter(Gene == "STING") %>%
    select(Probe, Project, BetaValue) %>%
    pivot_wider(
      names_from = Project,
      values_from = BetaValue,
      values_fn = mean  # Use mean if there are duplicates
    ) %>%
    select(Probe, all_of(order_vec))
  
  return(list(CGAS = cgas_table, STING = sting_table))
}

# Create the tables, passing the new column order
tables <- create_gene_tables(plot_data, project_col_order)

# Create results directory if it does not exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Save as an Excel file
output_file <- "D:/AcMet/3AcMet/epigenetic-analysis/Methylation/results/methylation_beta_values_ordered.xlsx"
write_xlsx(tables, output_file)

cat("\n[SUCCESS] Excel file successfully created:", output_file, "\n")