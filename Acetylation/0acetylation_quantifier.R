library(rtracklayer)
library(GenomicRanges)
library(writexl)
library(dplyr)
library(tidyr)

# Funzione principale per calcolare il segnale di acetilazione
calculate_acetylation_signal <- function(bigwig_path, chrom, start, end) {
  # Creazione dell'oggetto GRanges per la regione di interesse
  region_gr <- GRanges(seqnames = chrom,
                       ranges = IRanges(start = start, end = end))
  
  tryCatch({
    # Importazione del segnale dalla regione specificata
    signal_data <- import(bigwig_path, which = region_gr)
    
    # Se non ci sono dati nella regione
    if (length(signal_data) == 0) {
      return(0)  # Ritorna 0 invece di NA per evitare problemi nei confronti
    }
    
    # Calcolo del segnale medio normalizzato per la lunghezza della regione
    total_signal <- sum(signal_data$score)
    region_length <- end - start + 1
    normalized_signal <- total_signal / region_length
    
    return(normalized_signal)
  }, error = function(e) {
    message(paste("Errore nel processare il file:", bigwig_path))
    message(paste("Dettaglio errore:", e$message))
    return(NA)
  })
}

# Funzione per costruire la tabella dei risultati per un promotore
build_promoter_signal_table <- function(promoter_name, chrom, start, end, bigwig_files, cell_lines) {
  # Verifica che i file e le linee cellulari corrispondano
  if (length(bigwig_files) != length(cell_lines)) {
    stop("Il numero di file bigWig e linee cellulari deve essere uguale")
  }
  
  # Calcola il segnale per ogni file
  signals <- sapply(bigwig_files, function(bw_file) {
    calculate_acetylation_signal(bw_file, chrom, start, end)
  })
  
  # Crea il dataframe dei risultati
  results <- data.frame(
    Promoter = promoter_name,
    Cell_line = cell_lines,
    Mean_signal = signals,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

write_results_to_excel <- function(data, output_file) {
  
  # Estrai i nomi unici dei promotori
  promoters <- unique(data$Promoter)
  
  # Crea una lista per contenere i dataframe formattati (un foglio per promotore)
  sheet_list <- list()
  
  for (promoter in promoters) {
    
    # 1. Filtra i dati per il promotore corrente
    promoter_data <- data %>%
      filter(Promoter == promoter)
    
    # 2. Trasforma i dati nel formato 'wide' richiesto dallo script di plotting
    #    La colonna 'Probe' (che sarà il nome del promotore) e le colonne delle linee cellulari
    formatted_sheet <- promoter_data %>%
      select(Promoter, Cell_line, Mean_signal) %>%
      pivot_wider(
        names_from = Cell_line,
        values_from = Mean_signal
      ) %>%
      rename(Promoter = Promoter)
    
    # Inserisci il dataframe formattato nella lista con il nome del foglio = nome del promotore
    sheet_list[[promoter]] <- formatted_sheet
  }
  
  # Scrivi la lista di dataframe come fogli separati nel file Excel
  tryCatch({
    writexl::write_xlsx(sheet_list, path = output_file)
    cat("\n✓ Risultati salvati con successo in:", output_file, "\n")
  }, error = function(e) {
    cat("\n✗ Errore nel salvataggio del file Excel:", e$message, "\n")
  })
}

# Funzione per scrivere i risultati Mean_signal (Mean FC) in un singolo foglio wide
write_mean_fc_to_single_sheet_excel <- function(data, output_file) {
  
  # Trasforma i dati dal formato 'long' a 'wide'
  formatted_sheet <- data %>%
    select(Promoter, Cell_line, Mean_signal) %>%
    pivot_wider(
      names_from = Cell_line,
      values_from = Mean_signal
    ) %>%
    # Rinominare le colonne per corrispondere esattamente al tuo esempio
    rename(Promotore = Promoter) %>%
    # Aggiunge il suffisso (fold-change) ai nomi delle colonne delle linee cellulari
    rename_with(~ paste0(.x), .cols = -Promotore)
  
  # Scrive in un unico foglio chiamato "Mean_Fold_Change"
  sheet_list <- list(Mean_Fold_Change = formatted_sheet)
  
  tryCatch({
    writexl::write_xlsx(sheet_list, path = output_file)
    cat("\n✓ Tabella Fold-Change Media salvata con successo in:", output_file, "\n")
  }, error = function(e) {
    cat("\n✗ Errore nel salvataggio del file Excel:", e$message, "\n")
  })
}


# Esempio di utilizzo per CGAS e STING1
# Definisci i percorsi dei file BigWig (esempio)
bigwig_files <- c(
  "D:/AcMet/3AcMet/Epigenetic_Analysis/Acetylation/data/ENCFF681WFO.bigWig", #A549
  "D:/AcMet/3AcMet/Epigenetic_Analysis/Acetylation/data/ENCFF984WLE.bigWig"  #HCT116
)

#cell_lines <- c("A549", "HCT116", "LNCaP", "Caco-2")
cell_lines <- c("A549", "HCT116")

# Calcola per il promotore di CGAS
cgas_results <- build_promoter_signal_table(
  "cGAS", "chr6", 73450297, 73452797, 
  bigwig_files, cell_lines
)

# Calcola per il promotore di STING1
sting_results <- build_promoter_signal_table(
  "STING", "chr5", 139480935, 139483435, 
  bigwig_files, cell_lines
)

# Combina i risultati
final_results <- rbind(cgas_results, sting_results)
print(final_results)

# Aggiungi qui l'output di destinazione (usando un nome e percorso sensato per i risultati)
output_excel_file <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Acetylation/results/acetylation_signal_results.xlsx" 

# Verifica e crea la directory 'results' se necessario (assumendo una struttura simile al tuo plotting script)
output_dir <- dirname(output_excel_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Directory di output creata:", output_dir, "\n")
}

# Chiama la funzione di salvataggio
write_results_to_excel(final_results, output_excel_file)

output_excel_file_single <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Acetylation/results/acetylation_signal_results_single.xlsx"

write_mean_fc_to_single_sheet_excel(final_results, output_excel_file_single)
