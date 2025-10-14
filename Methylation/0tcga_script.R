# ANALISI METILAZIONE CGAS E STING IN DATASET TCGA - VERSIONE DEFINITIVA
setwd("D:/AcMet/3AcMet/Epigenetic_Analysis/Methylation/")

# CARICAMENTO PACCHETTI
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
  cat("✓ Pacchetti caricati correttamente\n")
}

# FUNZIONE PRINCIPALE OTTIMIZZATA
load_existing_methylation_data_final <- function() {
  cat("=== CARICAMENTO DATI METILAZIONE - VERSIONE FINALE ===\n")
  
  # 1. Definizione sonde target
  cgas_probes <- c("cg00996469", "cg09931909", "cg021283926", "cg16353624", 
                   "cg08905652", "cg16390139", "cg27279904", "cg00463577")
  sting_probes <- c("cg04232128", "cg03317505", "cg04560810", "cg14655316", 
                    "cg16532438")
  target_probes <- unique(c(cgas_probes, sting_probes))
  cat("Sonde target cercate:", length(target_probes), "\n")
  
  # 2. Trova i file .sesame.level3betas.txt
  data_dir <- "GDCdata"
  file_paths <- list.files(
    path = data_dir, 
    pattern = "\\.sesame\\.level3betas\\.txt$", 
    recursive = TRUE, 
    full.names = TRUE
  )
  
  cat("Trovati", length(file_paths), "file .sesame.level3betas.txt\n")
  
  if (length(file_paths) == 0) {
    stop("Nessun file '.sesame.level3betas.txt' trovato")
  }
  
  # 3. Funzione per estrarre TCGA barcode dal percorso
  extract_barcode <- function(path) {
    cat("Analizzo percorso:", path, "\n")
    
    # Metodo 1: Cerca nel percorso il pattern TCGA-XX-XXXX-XX
    barcode_match <- regmatches(path, regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]", path))
    if (length(barcode_match) > 0) {
      cat("Barcode trovato con regex:", barcode_match[1], "\n")
      return(barcode_match[1])
    }
    
    # Metodo 2: Estrai dalle directory
    path_parts <- strsplit(path, "/")[[1]]
    cat("Parti del percorso:", paste(path_parts, collapse = " | "), "\n")
    
    # Cerca parti che contengono TCGA
    tcga_parts <- path_parts[grepl("TCGA", path_parts)]
    if (length(tcga_parts) > 0) {
      # Prendi l'ultima parte che contiene TCGA (più specifica)
      potential_barcode <- tcga_parts[length(tcga_parts)]
      cat("Potenziale barcode da directory:", potential_barcode, "\n")
      return(potential_barcode)
    }
    
    cat("Nessun barcode trovato per:", path, "\n")
    return(NA)
  }
  
  # 4. Leggi i file con debugging
  cat("Inizio lettura dei file...\n")
  beta_list <- list()
  sample_barcodes <- character()
  files_processed <- 0
  files_with_data <- 0
  
  for (i in seq_along(file_paths)) {
    if (i %% 50 == 0) cat("Processati", i, "file su", length(file_paths), "\n")
    
    # Estrai barcode
    barcode <- extract_barcode(file_paths[i])
    
    if (!is.na(barcode)) {
      # Leggi il file
      data <- tryCatch({
        read.delim(file_paths[i], header = FALSE, stringsAsFactors = FALSE, sep = "\t")
      }, error = function(e) {
        cat("❌ Errore lettura file", file_paths[i], "-", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(data) && ncol(data) == 2) {
        colnames(data) <- c("Probe", "BetaValue")
        
        # Verifica che contenga le sonde target
        has_target_probes <- any(target_probes %in% data$Probe)
        
        if (has_target_probes) {
          beta_vector <- data$BetaValue
          names(beta_vector) <- data$Probe
          beta_list[[barcode]] <- beta_vector
          sample_barcodes <- c(sample_barcodes, barcode)
          files_with_data <- files_with_data + 1
          cat("✓ File", i, "- Barcode:", barcode, "- Sonde target: SI\n")
        } else {
          cat("⚠️ File", i, "- Barcode:", barcode, "- Sonde target: NO\n")
        }
      }
      files_processed <- files_processed + 1
    }
  }
  
  cat("File processati:", files_processed, "\n")
  cat("File con dati validi:", files_with_data, "\n")
  cat("Barcode estratti:", length(sample_barcodes), "\n")
  
  if (length(beta_list) == 0) {
    stop("Nessun dato valido trovato nei file")
  }
  
  # 5. DEBUG: Controlla quali sonde sono presenti
  all_probes <- unique(unlist(lapply(beta_list, names)))
  cat("Tutte le sonde trovate nei file (prime 20):", head(all_probes, 20), "\n")
  cat("Numero totale di sonde uniche:", length(all_probes), "\n")
  
  available_probes <- target_probes[target_probes %in% all_probes]
  cat("Sonde target disponibili:", length(available_probes), "/", length(target_probes), "\n")
  cat("Sonde target trovate:", paste(available_probes, collapse = ", "), "\n")
  
  if (length(available_probes) == 0) {
    stop("Nessuna delle sonde target è presente nei dati")
  }
  
  # 6. Crea matrice beta values
  cat("Creazione matrice beta values...\n")
  beta_matrix <- matrix(
    NA, 
    nrow = length(available_probes), 
    ncol = length(beta_list),
    dimnames = list(available_probes, names(beta_list))
  )
  
  # Riempimento matrice
  for (barcode in names(beta_list)) {
    probe_data <- beta_list[[barcode]]
    # Usa match per allineare correttamente le sonde
    match_idx <- match(available_probes, names(probe_data))
    beta_matrix[available_probes, barcode] <- probe_data[match_idx]
  }
  
  cat("Matrice creata - Dimensioni:", dim(beta_matrix), "\n")
  
  # 7. Crea metadati
  cat("Creazione metadati...\n")
  metadata <- data.frame(
    submitter_id = names(beta_list),
    project_id = substr(names(beta_list), 1, 9),  # Es: "TCGA-COAD"
    sample_type = "Primary Tumor",
    stringsAsFactors = FALSE
  )
  
  cat("Metadati creati per", nrow(metadata), "campioni\n")
  
  return(list(
    beta_values = beta_matrix,
    metadata = metadata,
    available_probes = available_probes,
    cgas_probes = cgas_probes,
    sting_probes = sting_probes
  ))
}

# SALVA DATI INTERMEDI
save_intermediate_data <- function(data, step_name) {
  filename <- paste0("intermediate_", step_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
  save(data, file = filename)
  cat("✓ Dati salvati in", filename, "\n")
  return(filename)
}

# CARICA DATI SALVATI
load_saved_data <- function(filename) {
  if (file.exists(filename)) {
    load(filename)
    cat("✓ Dati caricati da", filename, "\n")
    return(data)
  } else {
    stop("File non trovato:", filename)
  }
}

# ANALISI DESCRITTIVA CON CONTROLLI
perform_descriptive_analysis <- function(data) {
  cat("=== ANALISI DESCRITTIVA ===\n")
  
  # Controlli preliminari
  cat("Controllo dimensioni dati:\n")
  cat("  Beta values:", dim(data$beta_values), "\n")
  cat("  Metadati:", nrow(data$metadata), "\n")
  cat("  Sonde disponibili:", length(data$available_probes), "\n")
  
  # Creazione tabella riassuntiva
  summary_table <- data.frame(
    Sample = colnames(data$beta_values),
    Project = data$metadata$project_id,
    stringsAsFactors = FALSE
  )
  
  # Aggiungi i beta values per ogni sonda
  for(probe in data$available_probes) {
    beta_vals <- as.numeric(data$beta_values[probe, ])
    cat("Sonda", probe, "- Valori non-NA:", sum(!is.na(beta_vals)), "/", length(beta_vals), "\n")
    summary_table[[probe]] <- beta_vals
  }
  
  # Calcolo statistiche descrittive
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
  
  print(stats)
  
  return(list(summary_table = summary_table, statistics = stats))
}

# CREAZIONE VISUALIZZAZIONI CON CONTROLLI
create_visualizations <- function(data, descriptive_results) {
  cat("=== CREAZIONE VISUALIZZAZIONI ===\n")
  
  # Prepara dati per plotting
  plot_data <- descriptive_results$summary_table %>%
    pivot_longer(
      cols = -c(Sample, Project),
      names_to = "Probe",
      values_to = "BetaValue"
    ) %>%
    mutate(
      Gene = ifelse(Probe %in% data$cgas_probes, "CGAS", "STING")
    )
  
  cat("Dati per plotting:", nrow(plot_data), "righe\n")
  cat("Valori non-NA:", sum(!is.na(plot_data$BetaValue)), "\n")
  
  # Boxplot per gene e progetto
  p1 <- ggplot(plot_data, aes(x = Project, y = BetaValue, fill = Gene)) +
    geom_boxplot(alpha = 0.8, na.rm = TRUE) +
    labs(title = "Metilazione CGAS e STING in Tumori Primari TCGA",
         subtitle = paste("Basato su", length(data$available_probes), "sonde specifiche"),
         y = "Valore Beta di Metilazione", 
         x = "Tipo Tumorale") +
    theme_minimal() +
    theme(legend.position = "top")
  
  print(p1)
  
  # Salva il plot
  ggsave("results/metilazione_cgas_sting_boxplot.png", p1, width = 10, height = 6)
  
  # Boxplot separati per sonda se non troppe
  if (length(data$available_probes) <= 15) {
    p2 <- ggplot(plot_data, aes(x = Project, y = BetaValue, fill = Project)) +
      geom_boxplot(na.rm = TRUE) +
      facet_wrap(~ Probe + Gene, scales = "free_y", ncol = 3) +
      labs(title = "Metilazione per Sonda Specifica",
           y = "Valore Beta di Metilazione") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p2)
    ggsave("results/metilazione_per_sonde.png", p2, width = 12, height = 8)
  }
  
  return(plot_data)
}

# ESPORTAZIONE RISULTATI
export_results <- function(descriptive_results, plot_data, data) {
  cat("=== ESPORTAZIONE RISULTATI ===\n")
  
  # Crea directory se non esistono
  if (!dir.exists("results")) {
    dir.create("results")
  }
  
  # File paths
  beta_file <- "results/methylation_beta_values_target_probes.csv"
  stats_file <- "results/methylation_statistics_target_probes.csv"
  plot_data_file <- "results/plot_data.csv"
  probe_info_file <- "results/probe_information.csv"
  
  # Esporta dati
  write.csv(descriptive_results$summary_table, beta_file, row.names = FALSE)
  write.csv(descriptive_results$statistics, stats_file, row.names = FALSE)
  write.csv(plot_data, plot_data_file, row.names = FALSE)
  
  # Crea e esporta informazioni sulle sonde
  probe_info <- data.frame(
    Probe = data$available_probes,
    Gene = ifelse(data$available_probes %in% data$cgas_probes, "CGAS", "STING"),
    Status = "Presente nei dati",
    stringsAsFactors = FALSE
  )
  write.csv(probe_info, probe_info_file, row.names = FALSE)
  
  cat("✓ File esportati nella cartella 'results/':\n")
  cat("  -", beta_file, "\n")
  cat("  -", stats_file, "\n")
  cat("  -", plot_data_file, "\n")
  cat("  -", probe_info_file, "\n")
}

# FUNZIONE PRINCIPALE CON GESTIONE ERRORI
analyze_methylation_target_probes <- function(use_saved = FALSE, saved_file = NULL, force_download = FALSE) {
  cat("=== ANALISI METILAZIONE CGAS E STING - VERSIONE DEFINITIVA ===\n\n")
  
  tryCatch({
    # 1. Carica pacchetti
    load_required_packages()
    
    # 2. Gestione download (SOLO se force_download = TRUE)
    if (force_download) {
      cat("=== DOWNLOAD NUOVI DATI ===\n")
      query <- GDCquery(
        project = c("TCGA-COAD", "TCGA-LUAD", "TCGA-PRAD"),
        data.category = "DNA Methylation",
        data.type = "Methylation Beta Value",
        platform = "Illumina Human Methylation 450",
        sample.type = "Primary Tumor"
      )
      GDCdownload(query)
      cat("✓ Download completato\n")
    }
    
    # 3. Carica dati
    if (use_saved && !is.null(saved_file)) {
      cat("Caricamento dati salvati...\n")
      methylation_data <- load_saved_data(saved_file)
    } else {
      cat("Caricamento dati dai file locali...\n")
      methylation_data <- load_existing_methylation_data_final()
      
      # Salva i dati intermedi
      saved_file <- save_intermediate_data(methylation_data, "methylation_data")
      cat("Dati primari salvati in:", saved_file, "\n")
    }
    
    # 4. Analisi descrittiva
    descriptive_results <- perform_descriptive_analysis(methylation_data)
    save_intermediate_data(descriptive_results, "descriptive_analysis")
    
    # 5. Visualizzazioni
    plot_data <- create_visualizations(methylation_data, descriptive_results)
    save_intermediate_data(plot_data, "plot_data")
    
    # 6. Esportazione
    export_results(descriptive_results, plot_data, methylation_data)
    
    cat("\n🎉 === ANALISI COMPLETATA CON SUCCESSO === 🎉\n")
    cat("Dati salvati in:", saved_file, "\n")
    cat("Puoi riutilizzare i dati con: analyze_methylation_target_probes(use_saved = TRUE, saved_file = '", saved_file, "')\n")
    
    return(list(
      methylation_data = methylation_data,
      descriptive_results = descriptive_results,
      plot_data = plot_data,
      saved_file = saved_file
    ))
    
  }, error = function(e) {
    cat("❌ ERRORE durante l'analisi:", e$message, "\n")
    cat("Traccia dell'errore:\n")
    print(traceback())
    return(NULL)
  })
}

# ISTRUZIONI USO:
# - Per usare dati già scaricati: analyze_methylation_target_probes(use_saved = FALSE, force_download = FALSE)
# - Per nuovo download: analyze_methylation_target_probes(use_saved = FALSE, force_download = TRUE)
# - Per usare dati salvati: analyze_methylation_target_probes(use_saved = TRUE, saved_file = "tuo_file.RData")

# SPIEGAZIONE USO
cat("
ISTRUZIONI USO:

1. PRIMA ESECUZIONE:
   results <- analyze_methylation_target_probes()

2. ESECUZIONI SUCCESSIVE (con dati già processati):
   results <- analyze_methylation_target_probes(use_saved = TRUE, saved_file = 'intermediate_methylation_data_YYYYMMDD_HHMMSS.RData')

3. DEBUGGING:
   - I dati vengono salvati automaticamente dopo ogni step
   - Usa i file intermediate_* per debuggare
   - Controlla la cartella 'results/' per i file esportati
\n")

# ESEGUI L'ANALISI (SCEGLI UN'OPZIONE):

# Opzione 1: Esegui analisi completa (raccomandato per prima volta)
cat("Inizio analisi completa...\n")
results <- analyze_methylation_target_probes(use_saved = FALSE, force_download = FALSE)

# Opzione 2: Se hai già dei dati salvati, usa:
# results <- analyze_methylation_target_probes(use_saved = TRUE, saved_file = "intermediate_methylation_data_YYYYMMDD_HHMMSS.RData")

# Opzione 3: Solo caricamento dati (per debugging)
# methylation_data <- load_existing_methylation_data_final()
# save_intermediate_data(methylation_data, "debug_methylation_data")