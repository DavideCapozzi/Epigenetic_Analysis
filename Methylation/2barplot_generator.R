library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to process and visualize methylation data from each Excel sheet
process_methylation_sheet <- function(file_path, sheet_name, output_dir) {
  
  # Read the Excel sheet
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Ensure the first column is named "Probe"
  if (colnames(data)[1] != "Probe") {
    stop("The first column must be named 'Probe'")
  }
  
  # Extract the name of the Probe column
  probe_col <- colnames(data)[1]
  
  # Reshape data to long format for ggplot2
  data_long <- data %>%
    pivot_longer(
      cols = -all_of(probe_col),
      names_to = "Tumor_Type",
      values_to = "Beta_Value"
    ) %>%
    rename(Probe = all_of(probe_col)) %>%
    # Set Tumor_Type factor with the desired display order: COAD, PRAD, LUAD
    mutate(Tumor_Type = factor(Tumor_Type, levels = c("TCGA-COAD", "TCGA-PRAD", "TCGA-LUAD"))) %>%
    # Set Probe as a factor to preserve the correct order on the X axis
    mutate(Probe = factor(Probe, levels = unique(Probe)))
  
  # Labels mapping from TCGA codes to short tumor names
  tumor_labels <- c(
    "TCGA-COAD" = "COAD",
    "TCGA-LUAD" = "LUAD",
    "TCGA-PRAD" = "PRAD"
  )
  
  # Build the plot: Probe on X axis, faceted by Tumor Type
  p <- ggplot(data_long, aes(x = Probe, y = Beta_Value)) +
    # All bars filled with solid black
    geom_col(
      position = position_identity(),
      fill  = "black",   
      color = "black",   
      linewidth = 0.3,
      width = 0.75
    ) +
    
    # Y axis: beta-value range 0-1 with 0.1 breaks
    scale_y_continuous(
      name   = "β-value",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    
    # X axis: probe names
    scale_x_discrete(
      name = "Probe"
    ) +
    
    # Facet panels in the order COAD, PRAD, LUAD
    facet_wrap(
      ~ Tumor_Type,
      scales   = "free_x",
      ncol     = 3,
      labeller = as_labeller(tumor_labels)
    ) +
    
    # Plot title = sheet name
    ggtitle(sheet_name) +
    
    # Clean theme: no background grid lines
    theme_classic() +
    theme(
      # ALL TEXT SIZES INCREASED FOR READABILITY
      # X axis text: angled for readability, Arial 14
      axis.text.x  = element_text(angle = 45, hjust = 1,
                                  size = 14, family = "Arial"),
      # Y axis title: Arial 14 bold
      axis.title.y = element_text(size = 14, face = "bold", family = "Arial"),
      # Y axis tick labels: Arial 14
      axis.text.y  = element_text(size = 14, family = "Arial"),
      # X axis title: Arial 14 bold
      axis.title.x = element_text(size = 14, face = "bold", family = "Arial"),
      
      # Facet strip: no background box, bold Arial 14
      strip.background = element_blank(),
      strip.text       = element_text(size = 14, face = "bold", family = "Arial"),
      
      # No legend
      legend.position = "none",
      
      # Remove all panel grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Space between facet panels
      panel.spacing = unit(1.5, "cm"),
      
      # Plot title: INCREASED to 16 for better hierarchy
      plot.title = element_text(size = 16, face = "bold",
                                hjust = 0.5, family = "Arial")
    )
  
  # Build the output file path
  output_file <- file.path(output_dir, paste0(sheet_name, "_methylation.png"))
  
  # Save the plot to disk
  ggsave(output_file, plot = p, width = 12, height = 6, dpi = 300)
  
  cat("Plot saved:", output_file, "\n")
  
  return(p)
}

# MAIN - Set your file paths here before running
main <- function() {
  
  # ===== SET YOUR PARAMETERS HERE =====
  excel_file       <- "D:/AcMet/3AcMet/epigenetic-analysis/Methylation/results/methylation_beta_values_ordered.xlsx"
  output_directory <- "D:/AcMet/3AcMet/epigenetic-analysis/Methylation/imgs"
  
  # Check that the Excel file exists
  if (!file.exists(excel_file)) {
    stop("Excel file not found: ", excel_file)
  }
  
  # Create the output directory if it does not exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    cat("Output directory created:", output_directory, "\n")
  }
  
  # Read all sheet names from the workbook
  sheet_names <- excel_sheets(excel_file)
  cat("Sheets found:", paste(sheet_names, collapse = ", "), "\n\n")
  
  # Process each sheet
  for (sheet in sheet_names) {
    cat("Processing:", sheet, "\n")
    tryCatch(
      {
        process_methylation_sheet(excel_file, sheet, output_directory)
      },
      error = function(e) {
        cat("ERROR in sheet", sheet, ":", e$message, "\n")
      }
    )
  }
  
  cat("\nAll sheets processed.\n")
}

# Run the main function
main()