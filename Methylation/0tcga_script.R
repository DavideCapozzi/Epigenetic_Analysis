# METHYLATION ANALYSIS FOR CGAS AND STING IN TCGA DATASET - FINAL VERSION
setwd("D:/AcMet/3AcMet/epigenetic-analysis/Methylation/")

# LOAD REQUIRED PACKAGES
load_required_packages <- function() {
  required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", 
                         "ggplot2", "tidyr", "stringr")
  
  for (pkg in required_packages) {
    if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
      if(pkg %in% c("TCGAbiolinks", "SummarizedExperiment")) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
      library(pkg, character.only = TRUE)
    }
  }
  cat("[SUCCESS] Packages loaded successfully\n")
}

# OPTIMIZED MAIN FUNCTION
load_existing_methylation_data_final <- function() {
  cat("=== METHYLATION DATA LOADING - FINAL VERSION ===\n")
  
  # 1. Define target probes
  cgas_probes <- c("cg00996469", "cg09931909", "cg021283926", "cg16353624", 
                   "cg08905652", "cg16390139", "cg27279904", "cg00463577")
  sting_probes <- c("cg04232128", "cg03317505", "cg04560810", "cg14655316", 
                    "cg16532438")
  target_probes <- unique(c(cgas_probes, sting_probes))
  cat("[INFO] Target probes searched:", length(target_probes), "\n")
  
  # 2. Find .sesame.level3betas.txt files
  data_dir <- "GDCdata"
  file_paths <- list.files(
    path = data_dir, 
    pattern = "\\.sesame\\.level3betas\\.txt$", 
    recursive = TRUE, 
    full.names = TRUE
  )
  
  cat("[INFO] Found", length(file_paths), ".sesame.level3betas.txt files\n")
  
  if (length(file_paths) == 0) {
    stop("No '.sesame.level3betas.txt' files found")
  }
  
  # 3. Function to extract TCGA barcode from path
  extract_barcode <- function(path) {
    
    # Method 1: Search the path for the TCGA-XX-XXXX-XX pattern
    barcode_match <- regmatches(path, regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]", path))
    if (length(barcode_match) > 0) {
      return(barcode_match[1])
    }
    
    # Method 2: Extract from directories
    path_parts <- strsplit(path, "/")[[1]]
    
    # Search parts containing TCGA
    tcga_parts <- path_parts[grepl("TCGA", path_parts)]
    if (length(tcga_parts) > 0) {
      potential_barcode <- tcga_parts[length(tcga_parts)]
      return(potential_barcode)
    }
    
    cat("[WARN] No barcode found for:", path, "\n")
    return(NA)
  }
  
  # 4. Read files with debugging
  cat("[INFO] Starting file reading...\n")
  beta_list <- list()
  sample_barcodes <- character()
  files_processed <- 0
  files_with_data <- 0
  
  for (i in seq_along(file_paths)) {
    if (i %% 50 == 0) cat("[INFO] Processed", i, "files out of", length(file_paths), "\n")
    
    barcode <- extract_barcode(file_paths[i])
    
    if (!is.na(barcode)) {
      data <- tryCatch({
        read.delim(file_paths[i], header = FALSE, stringsAsFactors = FALSE, sep = "\t")
      }, error = function(e) {
        cat("[ERROR] Reading file", file_paths[i], "-", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(data) && ncol(data) == 2) {
        colnames(data) <- c("Probe", "BetaValue")
        
        has_target_probes <- any(target_probes %in% data$Probe)
        
        if (has_target_probes) {
          beta_vector <- data$BetaValue
          names(beta_vector) <- data$Probe
          beta_list[[barcode]] <- beta_vector
          sample_barcodes <- c(sample_barcodes, barcode)
          files_with_data <- files_with_data + 1
        }
      }
      files_processed <- files_processed + 1
    }
  }
  
  cat("[INFO] Files processed:", files_processed, "\n")
  cat("[INFO] Files with valid data:", files_with_data, "\n")
  cat("[INFO] Extracted barcodes:", length(sample_barcodes), "\n")
  
  if (length(beta_list) == 0) {
    stop("No valid data found in the files")
  }
  
  # 5. DEBUG: Check which probes are present
  all_probes <- unique(unlist(lapply(beta_list, names)))
  cat("[INFO] Total number of unique probes:", length(all_probes), "\n")
  
  available_probes <- target_probes[target_probes %in% all_probes]
  cat("[INFO] Target probes available:", length(available_probes), "/", length(target_probes), "\n")
  
  if (length(available_probes) == 0) {
    stop("None of the target probes are present in the data")
  }
  
  # 6. Create beta values matrix
  cat("[INFO] Creating beta values matrix...\n")
  beta_matrix <- matrix(
    NA, 
    nrow = length(available_probes), 
    ncol = length(beta_list),
    dimnames = list(available_probes, names(beta_list))
  )
  
  for (barcode in names(beta_list)) {
    probe_data <- beta_list[[barcode]]
    match_idx <- match(available_probes, names(probe_data))
    beta_matrix[available_probes, barcode] <- probe_data[match_idx]
  }
  
  cat("[INFO] Matrix created - Dimensions:", dim(beta_matrix), "\n")
  
  # 7. Create metadata
  cat("[INFO] Creating metadata...\n")
  metadata <- data.frame(
    submitter_id = names(beta_list),
    project_id = substr(names(beta_list), 1, 9), 
    sample_type = "Primary Tumor",
    stringsAsFactors = FALSE
  )
  
  cat("[INFO] Metadata created for", nrow(metadata), "samples\n")
  
  return(list(
    beta_values = beta_matrix,
    metadata = metadata,
    available_probes = available_probes,
    cgas_probes = cgas_probes,
    sting_probes = sting_probes
  ))
}

# SAVE INTERMEDIATE DATA
save_intermediate_data <- function(data, step_name) {
  filename <- paste0("intermediate_", step_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
  save(data, file = filename)
  cat("[SUCCESS] Data saved to", filename, "\n")
  return(filename)
}

# LOAD SAVED DATA
load_saved_data <- function(filename) {
  if (file.exists(filename)) {
    load(filename)
    cat("[SUCCESS] Data loaded from", filename, "\n")
    return(data)
  } else {
    stop("File not found:", filename)
  }
}

# DESCRIPTIVE ANALYSIS WITH CHECKS
perform_descriptive_analysis <- function(data) {
  cat("=== DESCRIPTIVE ANALYSIS ===\n")
  
  cat("[INFO] Checking data dimensions:\n")
  cat("  Beta values:", dim(data$beta_values), "\n")
  cat("  Metadata:", nrow(data$metadata), "\n")
  cat("  Available probes:", length(data$available_probes), "\n")
  
  summary_table <- data.frame(
    Sample = colnames(data$beta_values),
    Project = data$metadata$project_id,
    stringsAsFactors = FALSE
  )
  
  for(probe in data$available_probes) {
    beta_vals <- as.numeric(data$beta_values[probe, ])
    summary_table[[probe]] <- beta_vals
  }
  
  stats <- summary_table %>%
    group_by(Project) %>%
    summarise(
      n_samples = n(),
      across(where(is.numeric), 
             list(mean = ~mean(., na.rm = TRUE), 
                  sd = ~sd(., na.rm = TRUE),
                  median = ~median(., na.rm = TRUE),
                  na_count = ~sum(is.na(.))),
             .names = "{.col}_{.fn}")
    )
  
  return(list(summary_table = summary_table, statistics = stats))
}

# CREATE VISUALIZATIONS WITH CHECKS
create_visualizations <- function(data, descriptive_results) {
  cat("=== CREATING VISUALIZATIONS ===\n")
  
  plot_data <- descriptive_results$summary_table %>%
    pivot_longer(
      cols = -c(Sample, Project),
      names_to = "Probe",
      values_to = "BetaValue"
    ) %>%
    mutate(
      Gene = ifelse(Probe %in% data$cgas_probes, "CGAS", "STING")
    )
  
  cat("[INFO] Plotting data:", nrow(plot_data), "rows\n")
  cat("[INFO] Non-NA values:", sum(!is.na(plot_data$BetaValue)), "\n")
  
  p1 <- ggplot(plot_data, aes(x = Project, y = BetaValue, fill = Gene)) +
    geom_boxplot(alpha = 0.8, na.rm = TRUE) +
    labs(title = "CGAS and STING Methylation in TCGA Primary Tumors",
         subtitle = paste("Based on", length(data$available_probes), "specific probes"),
         y = "Methylation Beta Value", 
         x = "Tumor Type") +
    theme_minimal() +
    theme(legend.position = "top")
  
  if (!dir.exists("results")) {
    dir.create("results")
  }
  ggsave("results/methylation_cgas_sting_boxplot.png", p1, width = 10, height = 6)
  
  if (length(data$available_probes) <= 15) {
    p2 <- ggplot(plot_data, aes(x = Project, y = BetaValue, fill = Project)) +
      geom_boxplot(na.rm = TRUE) +
      facet_wrap(~ Probe + Gene, scales = "free_y", ncol = 3) +
      labs(title = "Methylation per Specific Probe",
           y = "Methylation Beta Value",
           x = "Tumor Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("results/methylation_per_probe.png", p2, width = 12, height = 8)
  }
  
  return(plot_data)
}

# EXPORT RESULTS
export_results <- function(descriptive_results, plot_data, data) {
  cat("=== EXPORTING RESULTS ===\n")
  
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  beta_file <- "results/methylation_beta_values_target_probes.csv"
  stats_file <- "results/methylation_statistics_target_probes.csv"
  plot_data_file <- "results/plot_data.csv"
  probe_info_file <- "results/probe_information.csv"
  
  write.csv(descriptive_results$summary_table, beta_file, row.names = FALSE)
  write.csv(descriptive_results$statistics, stats_file, row.names = FALSE)
  write.csv(plot_data, plot_data_file, row.names = FALSE)
  
  probe_info <- data.frame(
    Probe = data$available_probes,
    Gene = ifelse(data$available_probes %in% data$cgas_probes, "CGAS", "STING"),
    Status = "Present in data",
    stringsAsFactors = FALSE
  )
  write.csv(probe_info, probe_info_file, row.names = FALSE)
  
  cat("[SUCCESS] Files exported to the 'results/' folder.\n")
}

# MAIN FUNCTION WITH ERROR HANDLING
analyze_methylation_target_probes <- function(use_saved = FALSE, saved_file = NULL, force_download = FALSE) {
  cat("=== CGAS AND STING METHYLATION ANALYSIS - FINAL VERSION ===\n\n")
  
  tryCatch({
    load_required_packages()
    
    if (force_download) {
      cat("=== DOWNLOADING NEW DATA ===\n")
      query <- GDCquery(
        project = c("TCGA-COAD", "TCGA-LUAD", "TCGA-PRAD"),
        data.category = "DNA Methylation",
        data.type = "Methylation Beta Value",
        platform = "Illumina Human Methylation 450",
        sample.type = "Primary Tumor"
      )
      GDCdownload(query)
      cat("[SUCCESS] Download completed\n")
    }
    
    if (use_saved && !is.null(saved_file)) {
      cat("[INFO] Loading saved data...\n")
      methylation_data <- load_saved_data(saved_file)
    } else {
      cat("[INFO] Loading data from local files...\n")
      methylation_data <- load_existing_methylation_data_final()
      saved_file <- save_intermediate_data(methylation_data, "methylation_data")
      cat("[SUCCESS] Primary data saved to:", saved_file, "\n")
    }
    
    descriptive_results <- perform_descriptive_analysis(methylation_data)
    save_intermediate_data(descriptive_results, "descriptive_analysis")
    
    plot_data <- create_visualizations(methylation_data, descriptive_results)
    save_intermediate_data(plot_data, "plot_data")
    
    export_results(descriptive_results, plot_data, methylation_data)
    
    cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
    cat("[INFO] Data saved to:", saved_file, "\n")
    
    return(list(
      methylation_data = methylation_data,
      descriptive_results = descriptive_results,
      plot_data = plot_data,
      saved_file = saved_file
    ))
    
  }, error = function(e) {
    cat("[ERROR] During analysis:", e$message, "\n")
    return(NULL)
  })
}

# RUN THE ANALYSIS
cat("[INFO] Starting full analysis...\n")
results <- analyze_methylation_target_probes(use_saved = FALSE, force_download = FALSE)