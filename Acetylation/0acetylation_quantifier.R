library(rtracklayer)
library(GenomicRanges)
library(writexl)
library(dplyr)
library(tidyr)

# Main function to calculate the acetylation signal
calculate_acetylation_signal <- function(bigwig_path, chrom, start, end) {
  # Create the GRanges object for the region of interest
  region_gr <- GRanges(seqnames = chrom,
                       ranges = IRanges(start = start, end = end))
  
  tryCatch({
    # Import the signal from the specified region
    signal_data <- import(bigwig_path, which = region_gr)
    
    # If there is no data in the region
    if (length(signal_data) == 0) {
      return(0)  # Return 0 instead of NA to avoid comparison issues
    }
    
    # Calculate the mean signal normalized by the region length
    total_signal <- sum(signal_data$score)
    region_length <- end - start + 1
    normalized_signal <- total_signal / region_length
    
    return(normalized_signal)
  }, error = function(e) {
    message(paste("[ERROR] Processing file:", bigwig_path))
    message(paste("Error detail:", e$message))
    return(NA)
  })
}

# Function to build the results table for a specific promoter
build_promoter_signal_table <- function(promoter_name, chrom, start, end, bigwig_files, cell_lines) {
  # Verify that the number of files matches the number of cell lines
  if (length(bigwig_files) != length(cell_lines)) {
    stop("The number of bigWig files and cell lines must be equal")
  }
  
  # Calculate the signal for each file
  signals <- sapply(bigwig_files, function(bw_file) {
    calculate_acetylation_signal(bw_file, chrom, start, end)
  })
  
  # Create the results dataframe
  results <- data.frame(
    Promoter = promoter_name,
    Cell_line = cell_lines,
    Mean_signal = signals,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

write_results_to_excel <- function(data, output_file) {
  
  # Extract unique promoter names
  promoters <- unique(data$Promoter)
  
  # Create a list to hold formatted dataframes (one sheet per promoter)
  sheet_list <- list()
  
  for (promoter in promoters) {
    
    # 1. Filter data for the current promoter
    promoter_data <- data %>%
      filter(Promoter == promoter)
    
    # 2. Transform data into the 'wide' format required by the plotting script
    formatted_sheet <- promoter_data %>%
      select(Promoter, Cell_line, Mean_signal) %>%
      pivot_wider(
        names_from = Cell_line,
        values_from = Mean_signal
      ) 
    
    # Insert the formatted dataframe into the list
    sheet_list[[promoter]] <- formatted_sheet
  }
  
  # Write the list of dataframes as separate sheets in the Excel file
  tryCatch({
    writexl::write_xlsx(sheet_list, path = output_file)
    cat("\n[SUCCESS] Results successfully saved to:", output_file, "\n")
  }, error = function(e) {
    cat("\n[ERROR] Saving the Excel file:", e$message, "\n")
  })
}

# Function to write Mean_signal (Mean FC) results into a single wide sheet
write_mean_fc_to_single_sheet_excel <- function(data, output_file) {
  
  # Transform data from 'long' to 'wide' format
  formatted_sheet <- data %>%
    select(Promoter, Cell_line, Mean_signal) %>%
    pivot_wider(
      names_from = Cell_line,
      values_from = Mean_signal
    )
  
  # Write into a single sheet named "Mean_Fold_Change"
  sheet_list <- list(Mean_Fold_Change = formatted_sheet)
  
  tryCatch({
    writexl::write_xlsx(sheet_list, path = output_file)
    cat("\n[SUCCESS] Mean Fold-Change table successfully saved to:", output_file, "\n")
  }, error = function(e) {
    cat("\n[ERROR] Saving the Excel file:", e$message, "\n")
  })
}

# Example usage for cGAS and STING1
bigwig_files <- c(
  "D:/AcMet/3AcMet/epigenetic-analysis/Acetylation/data/ENCFF681WFO.bigWig", #A549
  "D:/AcMet/3AcMet/epigenetic-analysis/Acetylation/data/ENCFF984WLE.bigWig"  #HCT116
)

cell_lines <- c("A549", "HCT116")

# Calculate for the cGAS promoter
cgas_results <- build_promoter_signal_table(
  "cGAS", "chr6", 73450297, 73452797, 
  bigwig_files, cell_lines
)

# Calculate for the STING1 promoter
sting_results <- build_promoter_signal_table(
  "STING", "chr5", 139480935, 139483435, 
  bigwig_files, cell_lines
)

# Combine the results
final_results <- rbind(cgas_results, sting_results)
print(final_results)

# Define the destination output path
output_excel_file <- "D:/AcMet/3AcMet/epigenetic-analysis/Acetylation/results/acetylation_signal_results.xlsx" 

# Check and create the 'results' directory if necessary
output_dir <- dirname(output_excel_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("[INFO] Output directory created:", output_dir, "\n")
}

# Call the save functions
write_results_to_excel(final_results, output_excel_file)

output_excel_file_single <- "D:/AcMet/3AcMet/epigenetic-analysis/Acetylation/results/acetylation_signal_results_single.xlsx"

write_mean_fc_to_single_sheet_excel(final_results, output_excel_file_single)