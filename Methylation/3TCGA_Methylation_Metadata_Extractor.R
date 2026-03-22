# ==============================================================================
# SCRIPT: TCGA_Methylation_Metadata_Extractor_RStudio.R
# DESCRIZIONE: Estrae metadati essenziali dai file RData dell'analisi di metilazione
# USAGE: Esegui in RStudio - imposta rdata_path e chiama la funzione
# ==============================================================================

# CARICAMENTO PACCHETTI
load_required_packages <- function() {
  required_packages <- c("dplyr", "stringr")
  
  for (pkg in required_packages) {
    if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  cat("✅ Pacchetti caricati\n")
}

# FUNZIONE PRINCIPALE PER ESTRARRE METADATI
extract_methylation_metadata <- function(rdata_path = NULL) {
  cat("=== TCGA METHYLATION METADATA EXTRACTOR - RStudio ===\n\n")
  
  # 1. TROVA IL FILE RDATA
  if (is.null(rdata_path)) {
    cat("Nessun path specificato. Cerco l'ultimo file intermedio...\n")
    rdata_files <- list.files(
      pattern = "intermediate_methylation_data_.*\\.RData$",
      full.names = TRUE
    )
    
    if (length(rdata_files) == 0) {
      stop("❌ Nessun file RData trovato. Specifica il path manualmente.")
    }
    
    # Prendi l'ultimo file
    rdata_path <- sort(rdata_files, decreasing = TRUE)[1]
    cat("✅ Trovato file automaticamente:", rdata_path, "\n")
  } else {
    if (!file.exists(rdata_path)) {
      stop("❌ File non trovato: ", rdata_path)
    }
    cat("✅ File specificato:", rdata_path, "\n")
  }
  
  # 2. CARICA IL FILE RDATA
  cat("\n--- CARICAMENTO DATI ---\n")
  
  # Crea un ambiente temporaneo per il caricamento
  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  
  # Verifica quale oggetto è stato caricato
  loaded_objects <- ls(envir = temp_env)
  cat("Oggetti caricati:", paste(loaded_objects, collapse = ", "), "\n")
  
  # Cerca l'oggetto dei dati di metilazione (può avere nomi diversi)
  if ("methylation_data" %in% loaded_objects) {
    data <- temp_env$methylation_data
    cat("✅ Trovato oggetto 'methylation_data'\n")
  } else if ("data" %in% loaded_objects) {
    data <- temp_env$data
    cat("✅ Trovato oggetto 'data'\n")
  } else if (length(loaded_objects) == 1) {
    # Se c'è solo un oggetto, usalo
    data <- get(loaded_objects[1], envir = temp_env)
    cat("✅ Usato oggetto:", loaded_objects[1], "\n")
  } else {
    stop("❌ Non riesco a identificare l'oggetto con i dati di metilazione. Oggetti disponibili: ", 
         paste(loaded_objects, collapse = ", "))
  }
  
  # Verifica la struttura dei dati
  required_components <- c("beta_values", "metadata", "available_probes", "cgas_probes", "sting_probes")
  missing_components <- setdiff(required_components, names(data))
  
  if (length(missing_components) > 0) {
    stop("❌ Struttura dati incompleta. Componenti mancanti: ", paste(missing_components, collapse = ", "))
  }
  
  cat("✅ Struttura dati verificata\n")
  
  # 3. ESTRAZIONE METADATI ESSENZIALI
  cat("\n--- METADATI CAMPIONI ---\n")
  
  # Conteggio campioni per progetto
  sample_summary <- data$metadata %>%
    group_by(project_id) %>%
    summarise(
      n_campioni = n(),
      .groups = 'drop'
    )
  
  cat("CAMPIONI PER TUMORE:\n")
  print(sample_summary)
  cat("TOTALE CAMPIONI:", sum(sample_summary$n_campioni), "\n")
  
  # 4. INFORMAZIONI SONDE
  cat("\n--- INFORMAZIONI SONDE ---\n")
  
  cgas_found <- sum(data$cgas_probes %in% data$available_probes)
  sting_found <- sum(data$sting_probes %in% data$available_probes)
  cgas_missing <- setdiff(data$cgas_probes, data$available_probes)
  sting_missing <- setdiff(data$sting_probes, data$available_probes)
  
  cat("SONDE CGAS:", cgas_found, "/", length(data$cgas_probes), "trovate\n")
  cat("SONDE STING:", sting_found, "/", length(data$sting_probes), "trovate\n")
  cat("SONDE TOTALI ANALIZZATE:", length(data$available_probes), "\n")
  
  if (length(cgas_missing) > 0) {
    cat("SONDE CGAS MANCANTI:", paste(cgas_missing, collapse = ", "), "\n")
  }
  if (length(sting_missing) > 0) {
    cat("SONDE STING MANCANTI:", paste(sting_missing, collapse = ", "), "\n")
  }
  
  # 5. QUALITÀ DATI
  cat("\n--- QUALITÀ DATI ---\n")
  
  beta_matrix <- data$beta_values
  total_values <- length(beta_matrix)
  na_count <- sum(is.na(beta_matrix))
  na_percentage <- round(na_count / total_values * 100, 2)
  
  # Campioni con dati completi
  samples_with_data <- colSums(!is.na(beta_matrix))
  probes_per_sample <- nrow(beta_matrix)
  
  complete_samples <- sum(samples_with_data == probes_per_sample)
  good_samples <- sum(samples_with_data >= 0.9 * probes_per_sample)
  poor_samples <- sum(samples_with_data < 0.5 * probes_per_sample)
  
  cat("Dimensioni matrice beta:", dim(beta_matrix)[1], "sonde x", dim(beta_matrix)[2], "campioni\n")
  cat("Valori totali:", total_values, "\n")
  cat("Valori mancanti:", na_count, "(", na_percentage, "%)\n")
  cat("Campioni con dati completi (100% sonde):", complete_samples, "\n")
  cat("Campioni con dati buoni (>90% sonde):", good_samples, "\n")
  cat("Campioni con dati scarsi (<50% sonde):", poor_samples, "\n")
  
  # 6. BARCODE E CAMPIONI
  cat("\n--- BARCODE CAMPIONI ---\n")
  cat("Primi 10 barcode:\n")
  print(head(data$metadata$submitter_id, 10))
  cat("... e", length(data$metadata$submitter_id) - 10, "altri\n")
  
  # 7. CREAZIONE REPORT RIASSUNTIVO
  cat("\n", strrep("=", 50), "\n")
  cat("REPORT RIASSUNTIVO PER METODI\n")
  cat(strrep("=", 50), "\n")
  
  projects <- paste(unique(data$metadata$project_id), collapse = ", ")
  total_samples <- sum(sample_summary$n_campioni)
  
  cat("• FONTE DATI: TCGA (", projects, ")\n")
  cat("• PIATTAFORMA: Illumina Human Methylation 450K\n")
  cat("• PREPROCESSING: Dati normalizzati con pipeline SeSAMe del GDC\n")
  cat("• CAMPIONI: ")
  for (i in 1:nrow(sample_summary)) {
    cat(sample_summary$project_id[i], "=", sample_summary$n_campioni[i])
    if (i < nrow(sample_summary)) cat(", ")
  }
  cat(" (totale:", total_samples, "campioni tumorali primari)\n")
  cat("• SONDE: ", cgas_found, " cGAS + ", sting_found, " STING = ", 
      length(data$available_probes), " sonde totali\n", sep = "")
  cat("• QUALITÀ DATI: ", na_percentage, "% valori mancanti\n", sep = "")
  
  # 8. ESPORTAZIONE METADATI (OPZIONALE)
  cat("\n--- ESPORTAZIONE METADATI ---\n")
  
  output_dir <- "metadata_export"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    cat("✅ Creata directory:", output_dir, "\n")
  }
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # File completo metadati
  metadata_file <- file.path(output_dir, paste0("complete_metadata_", timestamp, ".csv"))
  write.csv(data$metadata, metadata_file, row.names = FALSE)
  cat("✅ Metadati completi esportati in:", metadata_file, "\n")
  
  # File riassuntivo
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
  
  cat("✅ Riassunto analisi esportato in:", summary_file, "\n")
  
  # 9. RESTITUISCI I DATI STRUTTURATI
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
  
  cat("\n✅ Estrazione completata!\n")
  return(invisible(results))
}

# ==============================================================================
# ISTRUZIONI PER L'USO IN RSTUDIO
# ==============================================================================

cat("
📋 ISTRUZIONI PER L'USO IN RSTUDIO:

1. CARICA I PACCHETTI:
   load_required_packages()

2. IMPOSTA IL PERCORSO DEL FILE RDATA:
   rdata_path <- 'tuo/path/al/file.RData'

   Esempi:
   rdata_path <- 'intermediate_methylation_data_20231201_143022.RData'
   rdata_path <- 'D:/Analisi/TCGA/methylation_data.RData'

3. ESEGUI L'ESTRAZIONE:
   result <- extract_methylation_metadata(rdata_path)

4. SE NON SPECIFICHI IL PATH, CERCHERA' AUTOMATICAMENTE:
   result <- extract_methylation_metadata()

5. ACCEDI AI RISULTATI:
   print(result$sample_summary)    # Campioni per tumore
   print(result$probe_summary)     # Sommario sonde
   print(result$data_quality)      # Qualità dati
   result$sample_barcodes          # Lista completa barcode

📊 OUTPUT:
   • Report dettagliato in console
   • File CSV con metadati completi
   • File di testo con riassunto
   • Oggetto R strutturato con tutti i dati

")

# ==============================================================================
# ESEMPIO DI ESECUZIONE (RIMUOVI IL COMMENTO PER PROVARE)
# ==============================================================================

# ESEMPIO 1: Con path specifico
rdata_path <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Methylation/intermediate_methylation_data_20251008_215333.RData"
result <- extract_methylation_metadata(rdata_path)

# ESEMPIO 2: Ricerca automatica
# result <- extract_methylation_metadata()

# ESEMPIO 3: Usa i risultati
# print(result$sample_summary)
# cat("Campioni totali:", result$data_quality$total_samples, "\n")

cat("🚀 PRONTO PER L'USO! Imposta rdata_path e chiama extract_methylation_metadata()\n")