library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Funzione per processare e visualizzare i dati da ogni foglio
process_methylation_sheet <- function(file_path, sheet_name, output_dir) {
  
  # Leggi il foglio Excel
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Verifica che la prima colonna sia Probe
  if (colnames(data)[1] != "Probe") {
    stop("La prima colonna deve essere 'Probe'")
  }
  
  # Estrai il nome della colonna Probe
  probe_col <- colnames(data)[1]
  
  # Trasforma in formato long per ggplot2
  data_long <- data %>%
    pivot_longer(
      cols = -all_of(probe_col),
      names_to = "Tumor_Type",
      values_to = "Beta_Value"
    ) %>%
    rename(Probe = all_of(probe_col)) %>% 
    # Converti Tumor_Type a factor per controllare l'ordine
    mutate(Tumor_Type = factor(Tumor_Type, levels = unique(Tumor_Type))) %>%
    # Aggiungi Probe come fattore per assicurare l'ordine corretto sull'asse X
    mutate(Probe = factor(Probe, levels = unique(Probe)))
  
  tumor_labels <- c(
    "TCGA-COAD" = "COAD", # Esempio di abbreviazione
    "TCGA-LUAD" = "LUAD",
    "TCGA-PRAD" = "PRAD"
  )
  
  # Prepara la mappatura dei colori (stessa logica di prima)
  tumor_levels <- levels(data_long$Tumor_Type)
  color_values <- c("#ECDD90", "#8BA345", "#6CADDE") # Ho mantenuto i tuoi colori originali
  
  color_mapping <- color_values[1:min(length(tumor_levels), length(color_values))]
  names(color_mapping) <- tumor_levels[1:min(length(tumor_levels), length(color_values))]
  
  # Crea il grafico: Probe sull'asse X, Facciata (Facet) per il Tipo di Tumore
  p <- ggplot(data_long, aes(x = Probe, y = Beta_Value)) + 
    # Barre adiacenti all'interno di ogni facciata, con colore per Tipo di Tumore
    geom_col(position = position_identity(), 
             color = "black", 
             linewidth = 0.3,
             width = 0.75, 
             aes(fill = Tumor_Type)) + # <-- Fill specificato qui
    
    scale_y_continuous(
      name = "β-value",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    scale_x_discrete(
      name = ""
    ) +
    # APPLICA LA SEPARAZIONE DEI GRUPPI TRAMITE FACET_WRAP
    facet_wrap(~ Tumor_Type, 
               scales = "free_x", # Assicura che l'asse X sia separato
               ncol = length(tumor_levels), # Disponi i pannelli su una sola riga
               labeller = as_labeller(tumor_labels)) +
    
    scale_fill_manual(
      values = color_mapping
    ) +
    ggtitle(sheet_name) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title.y = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      
      # Elementi che controllano il Facet (Etichetta e Sfondo)
      strip.background = element_blank(), # Rimuove lo sfondo del titolo della facciata
      strip.text = element_text(size = 11, face = "bold"), # Stile per il titolo del tumore
      
      # Legenda
      legend.position = "none", # Rimuovi la legenda se il titolo del facet è sufficiente
      
      #panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(1.5, "cm"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # Crea il nome del file di output
  output_file <- file.path(output_dir, paste0(sheet_name, "_methylation.png"))
  
  # Salva il grafico
  ggsave(output_file, plot = p, width = 12, height = 6, dpi = 300)
  
  cat("✓ Grafico salvato:", output_file, "\n")
  
  return(p)
}

# MAIN - Specifica qui i tuoi parametri
main <- function() {
  
  # ===== SPECIFICA I TUOI PARAMETRI QUI =====
  excel_file <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Methylation/results/methylation_beta_values_ordered.xlsx"  # Sostituisci con il percorso del tuo file
  output_directory <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Methylation/imgs"     # Sostituisci con la directory di output
  
  # Verifica che il file Excel esista
  if (!file.exists(excel_file)) {
    stop("Il file Excel non esiste: ", excel_file)
  }
  
  # Crea la directory di output se non esiste
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    cat("Directory di output creata:", output_directory, "\n")
  }
  
  # Leggi i nomi di tutti i fogli
  sheet_names <- excel_sheets(excel_file)
  cat("Fogli trovati:", paste(sheet_names, collapse = ", "), "\n\n")
  
  # Processa ogni foglio
  for (sheet in sheet_names) {
    cat("Processando:", sheet, "\n")
    tryCatch(
      {
        process_methylation_sheet(excel_file, sheet, output_directory)
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
