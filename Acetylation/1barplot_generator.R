library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 
library(scales) # Necessario per formattare i numeri nelle etichette

# Funzione per processare e visualizzare i dati da ogni foglio
process_acetylation_sheet <- function(file_path, sheet_name, output_dir) {
  
  # Leggi il foglio Excel
  data <- read_excel(file_path, sheet = sheet_name)
  
  # La prima colonna deve essere 'Promoter' (come concordato per l'Excel semplificato)
  if (colnames(data)[1] != "Promoter") {
    stop("La prima colonna del foglio Excel deve essere 'Promoter'.")
  }
  
  # Estrai il nome della colonna
  promoter_col <- colnames(data)[1] 
  
  # Trasforma in formato long per ggplot2
  data_long <- data %>%
    pivot_longer(
      cols = -all_of(promoter_col),
      names_to = "Cell_Line",
      values_to = "Mean_Signal" 
    ) %>%
    rename(Promoter = all_of(promoter_col)) %>%  
    # Converti a factor per controllare l'ordine
    mutate(Cell_Line = factor(Cell_Line, levels = unique(Cell_Line))) %>%
    mutate(Promoter = factor(Promoter, levels = unique(Promoter)))
  
  # ----------------------------------------------------
  # Colori
  # ----------------------------------------------------
  cell_levels <- levels(data_long$Cell_Line)
  color_values <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")  
  
  color_mapping <- color_values[1:min(length(cell_levels), length(color_values))]
  names(color_mapping) <- cell_levels[1:min(length(cell_levels), length(color_values))]
  
  # ----------------------------------------------------
  # CREAZIONE DEL GRAFICO BARPLOT MODIFICATO
  # ----------------------------------------------------
  p <- ggplot(data_long, aes(x = Promoter, y = Mean_Signal, fill = Cell_Line)) + 
    # Barre
    geom_col(position = position_dodge(width = 0.8), 
             color = "black", 
             linewidth = 0.3,
             width = 0.75) + 
    
    # Etichette dei valori sopra le barre
    geom_text(
      aes(label = format(round(Mean_Signal, 2), nsmall = 2)),
      position = position_dodge(width = 0.8), 
      vjust = -0.5,
      size = 3.5,
      angle = 0, # Manteniamo l'etichetta orizzontale se possibile
      hjust = 0.5 
    ) +
    
    # Facet per Linea Cellulare, etichette in basso, senza titolo
    facet_wrap(~ Cell_Line, 
               scales = "free_x", 
               ncol = length(cell_levels), 
               strip.position = "bottom") + 
    
    scale_y_continuous(
      name = "Mean Acetylation Signal (BigWig Score)",
      limits = c(0, max(data_long$Mean_Signal) * 1.2), # Aumentato leggermente lo spazio per le etichette orizzontali
      breaks = pretty_breaks(n = 10)
    ) +
    scale_x_discrete(
      name = "" # Etichetta dell'asse X principale vuota
    ) +
    
    scale_fill_manual(
      values = color_mapping
    ) +
    ggtitle(sheet_name) +
    theme_minimal() +
    theme(
      # RIMOZIONE ETICHETTE DEI PROMOTORI SOTTO LE BARRE
      axis.text.x = element_blank(), # Rimuove le etichette (nome del gene)
      axis.ticks.x = element_blank(), # Rimuove i piccoli segni di spunta
      
      axis.title.y = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      
      # Rimuove lo sfondo dei facet e il titolo (ma lascia l'etichetta in fondo)
      strip.background = element_blank(), 
      strip.text = element_text(size = 11, face = "bold"), 
      
      legend.position = "none",
      
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(0.5, "cm"), 
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      
      # L'etichetta del Facet (Cell_Line) è ora l'unica informazione sull'asse X
      axis.title.x = element_blank()
    )
  
  # Crea il nome del file di output
  output_file <- file.path(output_dir, paste0(sheet_name, "_acetylation.png")) 
  
  # Salva il grafico
  ggsave(output_file, plot = p, width = 10, height = 7, dpi = 300)
  
  cat("✓ Grafico salvato:", output_file, "\n")
  
  return(p)
}

# MAIN - Specifica qui i tuoi parametri
main <- function() {
  
  # ===== SPECIFICA I TUOI PARAMETRI QUI =====
  excel_file <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Acetylation/results/acetylation_signal_results.xlsx" 
  output_directory <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Acetylation/imgs"    
  
  # Verifica e crea le directory
  if (!file.exists(excel_file)) {
    stop("Il file Excel non esiste: ", excel_file)
  }
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    cat("Directory di output creata:", output_directory, "\n")
  }
  
  # Leggi i nomi di tutti i fogli
  sheet_names <- readxl::excel_sheets(excel_file)
  cat("Fogli trovati:", paste(sheet_names, collapse = ", "), "\n\n")
  
  # Processa ogni foglio
  for (sheet in sheet_names) {
    cat("Processando:", sheet, "\n")
    tryCatch(
      {
        process_acetylation_sheet(excel_file, sheet, output_directory)
      },
      error = function(e) {
        cat("✗ Errore nel foglio", sheet, ":", e$message, "\n")
      }
    )
  }
  
  cat("\nProcessamento completato!\n")
}

# Esegui il main
main()