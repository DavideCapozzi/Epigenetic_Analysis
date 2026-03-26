# ==============================================================================
# SCRIPT: TCGA_Methylation_Metadata_Extractor_RStudio.R
# DESCRIPTION: Extracts essential metadata from RData files of methylation analysis
# USAGE: Run in RStudio - set rdata_path and call the function
# ==============================================================================

setwd("D:/AcMet/3AcMet/epigenetic-analysis/Methylation/")

# LOAD REQUIRED PACKAGES
load_required_packages <- function() {
  required_packages <- c("dplyr", "stringr")
  
  for (pkg in required_packages) {
    if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  cat("[SUCCESS] Packages loaded\n")
}

# MAIN FUNCTION TO EXTRACT METADATA
extract_methylation_metadata <- function(rdata_path = NULL) {
  cat("=== TCGA METHYLATION METADATA EXTRACTOR - RStudio ===\n\n")
  
  # 1. FIND THE RDATA FILE
  if (is.null(rdata_path)) {
    cat("[INFO] No path specified. Searching for the latest intermediate file...\n")
    rdata_files <- list.files(
      pattern = "intermediate_methylation_data_.*\\.RData$",
      full.names = TRUE
    )
    
    if (length(rdata_files) == 0) {
      stop("[ERROR] No RData file found. Please specify the path manually.")
    }
    
    # Get the latest file
    rdata_path <- sort(rdata_files, decreasing = TRUE)[1]
    cat("[SUCCESS] File automatically found:", rdata_path, "\n")
  } else {
    if (!file.exists(rdata_path)) {
      stop("[ERROR] File not found: ", rdata_path)
    }
    cat("[SUCCESS] Specified file:", rdata_path, "\n")
  }
  
  # 2. LOAD THE RDATA FILE
  cat("\n--- DATA LOADING ---\n")
  
  # Create a temporary environment for loading
  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  
  # Check which object was loaded
  loaded_objects <- ls(envir = temp_env)
  cat("[INFO] Loaded objects:", paste(loaded_objects, collapse = ", "), "\n")
  
  # Look for the methylation data object
  if ("methylation_data" %in% loaded_objects) {
    data <- temp_env$methylation_data
    cat("[SUCCESS] Found 'methylation_data' object\n")
  } else if ("data" %in% loaded_objects) {
    data <- temp_env$data
    cat("[SUCCESS] Found 'data' object\n")
  } else if (length(loaded_objects) == 1) {
    # If there's only one object, use it
    data <- get(loaded_objects[1], envir = temp_env)
    cat("[SUCCESS] Used object:", loaded_objects[1], "\n")
  } else {
    stop("[ERROR] Cannot identify the methylation data object. Available objects: ", 
         paste(loaded_objects, collapse = ", "))
  }
  
  # Verify the data structure
  required_components <- c("beta_values", "metadata", "available_probes", "cgas_probes", "sting_probes")
  missing_components <- setdiff(required_components, names(data))
  
  if (length(missing_components) > 0) {
    stop("[ERROR] Incomplete data structure. Missing components: ", paste(missing_components, collapse = ", "))
  }
  
  cat("[SUCCESS] Data structure verified\n")
  
  # 3. ESSENTIAL METADATA EXTRACTION
  cat("\n--- SAMPLE METADATA ---\n")
  
  # Sample count per project
  sample_summary <- data$metadata %>%
    group_by(project_id) %>%
    summarise(
      n_samples = n(),
      .groups = 'drop'
    )
  
  cat("SAMPLES PER TUMOR:\n")
  print(sample_summary)
  cat("TOTAL SAMPLES:", sum(sample_summary$n_samples), "\n")
  
  # 4. PROBES INFORMATION
  cat("\n--- PROBES INFORMATION ---\n")
  
  cgas_found <- sum(data$cgas_probes %in% data$available_probes)
  sting_found <- sum(data$sting_probes %in% data$available_probes)
  cgas_missing <- setdiff(data$cgas_probes, data$available_probes)
  sting_missing <- setdiff(data$sting_probes, data$available_probes)
  
  cat("CGAS PROBES:", cgas_found, "/", length(data$cgas_probes), "found\n")
  cat("STING PROBES:", sting_found, "/", length(data$sting_probes), "found\n")
  cat("TOTAL PROBES ANALYZED:", length(data$available_probes), "\n")
  
  if (length(cgas_missing) > 0) {
    cat("MISSING CGAS PROBES:", paste(cgas_missing, collapse = ", "), "\n")
  }
  if (length(sting_missing) > 0) {
    cat("MISSING STING PROBES:", paste(sting_missing, collapse = ", "), "\n")
  }
  
  # 5. DATA QUALITY
  cat("\n--- DATA QUALITY ---\n")
  
  beta_matrix <- data$beta_values
  total_values <- length(beta_matrix)
  na_count <- sum(is.na(beta_matrix))
  na_percentage <- round(na_count / total_values * 100, 2)
  
  samples_with_data <- colSums(!is.na(beta_matrix))
  probes_per_sample <- nrow(beta_matrix)
  
  complete_samples <- sum(samples_with_data == probes_per_sample)
  good_samples <- sum(samples_with_data >= 0.9 * probes_per_sample)
  poor_samples <- sum(samples_with_data < 0.5 * probes_per_sample)
  
  cat("Beta matrix dimensions:", dim(beta_matrix)[1], "probes x", dim(beta_matrix)[2], "samples\n")
  cat("Total values:", total_values, "\n")
  cat("Missing values:", na_count, "(", na_percentage, "%)\n")
  cat("Samples with complete data (100% probes):", complete_samples, "\n")
  cat("Samples with good data (>90% probes):", good_samples, "\n")
  cat("Samples with poor data (<50% probes):", poor_samples, "\n")
  
  # 6. SUMMARY REPORT CREATION
  cat("\n", strrep("=", 50), "\n")
  cat("METHODS SUMMARY REPORT\n")
  cat(strrep("=", 50), "\n")
  
  projects <- paste(unique(data$metadata$project_id), collapse = ", ")
  total_samples <- sum(sample_summary$n_samples)
  
  cat("* DATA SOURCE: TCGA (", projects, ")\n")
  cat("* PLATFORM: Illumina Human Methylation 450K\n")
  cat("* PREPROCESSING: Data normalized with GDC SeSAMe pipeline\n")
  cat("* SAMPLES: ")
  for (i in 1:nrow(sample_summary)) {
    cat(sample_summary$project_id[i], "=", sample_summary$n_samples[i])
    if (i < nrow(sample_summary)) cat(", ")
  }
  cat(" (total:", total_samples, "primary tumor samples)\n")
  cat("* PROBES: ", cgas_found, " cGAS + ", sting_found, " STING = ", 
      length(data$available_probes), " total probes\n", sep = "")
  cat("* DATA QUALITY: ", na_percentage, "% missing values\n", sep = "")
  
  # 7. METADATA EXPORT (OPTIONAL)
  cat("\n--- METADATA EXPORT ---\n")
  
  output_dir <- "metadata_export"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    cat("[INFO] Created directory:", output_dir, "\n")
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Complete metadata file
  metadata_file <- file.path(output_dir, paste0("complete_metadata_", timestamp, ".csv"))
  write.csv(data$metadata, metadata_file, row.names = FALSE)
  cat("[SUCCESS] Complete metadata exported to:", metadata_file, "\n")
  
  # Summary file
  summary_file <- file.path(output_dir, paste0("analysis_summary_", timestamp, ".txt"))
  sink(summary_file)
  cat("ANALYSIS SUMMARY - TCGA Methylation Data\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Source file:", rdata_path, "\n")
  cat("=========================================\n\n")
  
  cat("SAMPLE SUMMARY:\n")
  print(sample_summary)
  cat("\nTotal samples:", total_samples, "\n\n")
  
  cat("PROBES SUMMARY:\n")
  cat("cGAS probes found:", cgas_found, "/", length(data$cgas_probes), "\n")
  cat("STING probes found:", sting_found, "/", length(data$sting_probes), "\n")
  cat("Total probes analyzed:", length(data$available_probes), "\n")
  cat("Missing cGAS probes:", if(length(cgas_missing)>0) paste(cgas_missing, collapse=", ") else "None", "\n")
  cat("Missing STING probes:", if(length(sting_missing)>0) paste(sting_missing, collapse=", ") else "None", "\n\n")
  
  cat("DATA QUALITY:\n")
  cat("Missing values:", na_count, "/", total_values, "(", na_percentage, "%)\n")
  cat("Samples with complete data:", complete_samples, "\n")
  cat("Samples with good data (>90%):", good_samples, "\n")
  cat("Samples with poor data (<50%):", poor_samples, "\n")
  sink()
  
  cat("[SUCCESS] Analysis summary exported to:", summary_file, "\n")
  
  # 8. RETURN STRUCTURED DATA
  results <- list(
    sample_summary = sample_summary,
    probe_summary = list(
      cgas_found = cgas_found,
      cgas_total = length(data$cgas_probes),
      sting_found = sting_found, 
      sting_total = length(data$sting_probes),
      available_probes = data$available_probes,
      missing_cgas = cgas_missing,
      missing_sting = sting_missing
    ),
    data_quality = list(
      na_percentage = na_percentage,
      total_samples = total_samples,
      total_probes = length(data$available_probes),
      complete_samples = complete_samples,
      good_samples = good_samples,
      poor_samples = poor_samples
    ),
    file_used = rdata_path,
    sample_barcodes = data$metadata$submitter_id
  )
  
  cat("\n[SUCCESS] Extraction completed!\n")
  return(invisible(results))
}

# ==============================================================================
# EXECUTION EXAMPLE
# ==============================================================================

rdata_path <- "D:/AcMet/3AcMet/epigenetic-analysis/Methylation/intermediate_methylation_data_20251008_215333.RData"
result <- extract_methylation_metadata(rdata_path)