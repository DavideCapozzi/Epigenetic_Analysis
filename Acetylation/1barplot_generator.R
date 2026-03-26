library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales) 

# Function to process and visualize acetylation data from each Excel sheet
process_acetylation_sheet <- function(file_path, sheet_name, output_dir) {
  
  # Read the Excel sheet
  data <- read_excel(file_path, sheet = sheet_name)
  
  # The first column must be named 'Promoter'
  if (colnames(data)[1] != "Promoter") {
    stop("The first column of the Excel sheet must be named 'Promoter'.")
  }
  
  # Extract the name of the first column
  promoter_col <- colnames(data)[1]
  
  # Reshape data to long format for ggplot2
  data_long <- data %>%
    pivot_longer(
      cols      = -all_of(promoter_col),
      names_to  = "Cell_Line",
      values_to = "Mean_Signal"
    ) %>%
    rename(Promoter = all_of(promoter_col)) %>%
    # Convert to factors to control display order: HCT116 first, then A549
    mutate(Cell_Line = factor(Cell_Line, levels = c("HCT116", "A549"))) %>%
    mutate(Promoter  = factor(Promoter,  levels = unique(Promoter)))
  
  # Number of cell lines (used to set facet columns)
  cell_levels <- levels(data_long$Cell_Line)
  
  # Build the bar plot
  p <- ggplot(data_long, aes(x = Cell_Line, y = Mean_Signal)) +
    # All bars filled with solid black; no fill aesthetic mapping
    geom_col(
      position  = position_dodge(width = 0.8),
      fill      = "black",   
      color     = "black",   
      linewidth = 0.3,
      width     = 0.75
    ) +
    
    # Value labels above bars
    geom_text(
      aes(label = format(round(Mean_Signal, 2), nsmall = 2)),
      position = position_dodge(width = 0.8),
      vjust    = -0.5,
      size     = 6,       
      angle    = 0,
      hjust    = 0.5,
      family   = "Arial"  
    ) +
    
    # Facet by Cell Line
    facet_wrap(
      ~ Cell_Line,
      scales         = "free_x",
      ncol           = length(cell_levels),
      strip.position = "bottom"
    ) +
    
    # Y axis: auto range with some headroom for labels
    scale_y_continuous(
      name   = "Mean Acetylation Signal (BigWig Score)",
      limits = c(0, max(data_long$Mean_Signal) * 1.2),
      breaks = pretty_breaks(n = 10)
    ) +
    
    # X axis: title left blank (cell line name shown in facet strip below)
    scale_x_discrete(
      name = ""
    ) +
    
    # Plot title = sheet name
    ggtitle(sheet_name) +
    
    # Start from a blank slate
    theme_void() +
    theme(
      # Draw only the left (Y) and bottom (X) axis lines
      axis.line.y       = element_line(color = "black", linewidth = 0.5),
      axis.line.x       = element_line(color = "black", linewidth = 0.5),
      
      # Y axis ticks
      axis.ticks.y      = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      
      # Y axis title and tick labels
      axis.title.y = element_text(size = 14, face = "bold", family = "Arial",
                                  angle = 90, vjust = 0.5, margin = margin(r = 6)),
      axis.text.y  = element_text(size = 14, family = "Arial", hjust = 1,
                                  margin = margin(r = 3)),
      
      # Hide X axis tick marks and text
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      
      # Facet strip
      strip.background = element_blank(),
      strip.placement  = "outside",
      strip.text       = element_text(size = 14, face = "bold", family = "Arial",
                                      margin = margin(t = 4, b = 15)),
      
      # No legend
      legend.position = "none",
      
      # No grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Space between facet panels
      panel.spacing = unit(0.5, "cm"),
      
      # Plot title
      plot.title = element_text(size = 16, face = "bold",
                                hjust = 0.5, family = "Arial",
                                margin = margin(t = 15, b = 8)),
      
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Build the output file path
  output_file <- file.path(output_dir, paste0(sheet_name, "_acetylation.png"))
  
  # Save the plot to disk
  ggsave(output_file, plot = p, width = 10, height = 7, dpi = 300)
  
  cat("[INFO] Plot saved:", output_file, "\n")
  
  return(p)
}

# MAIN - Set your file paths here before running
main <- function() {
  
  # ===== SET YOUR PARAMETERS HERE =====
  excel_file       <- "D:/AcMet/3AcMet/epigenetic-analysis/Acetylation/results/acetylation_signal_results.xlsx"
  output_directory <- "D:/AcMet/3AcMet/epigenetic-analysis/Acetylation/imgs"
  
  # Check that the Excel file exists
  if (!file.exists(excel_file)) {
    stop("Excel file not found: ", excel_file)
  }
  
  # Create the output directory if it does not exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    cat("[INFO] Output directory created:", output_directory, "\n")
  }
  
  # Read all sheet names from the workbook
  sheet_names <- readxl::excel_sheets(excel_file)
  cat("[INFO] Sheets found:", paste(sheet_names, collapse = ", "), "\n\n")
  
  # Process each sheet
  for (sheet in sheet_names) {
    cat("[INFO] Processing:", sheet, "\n")
    tryCatch(
      {
        process_acetylation_sheet(excel_file, sheet, output_directory)
      },
      error = function(e) {
        cat("[ERROR] In sheet", sheet, ":", e$message, "\n")
      }
    )
  }
  
  cat("\n[INFO] All sheets processed.\n")
}

# Run the main function
main()