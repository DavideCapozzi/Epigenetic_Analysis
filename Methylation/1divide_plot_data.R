# SCRIPT PER CREARE FILE EXCEL DA plot_data.csv
# Salva questo codice in un file separato, es. create_excel_tables.R

# Carica librerie necessarie
if (!require("writexl")) {
  install.packages("writexl")
  library(writexl)
}

if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require("tidyr")) {
  install.packages("tidyr")
  library(tidyr)
}

# Leggi i dati da plot_data.csv
# NOTA: Assicurati che il percorso sia corretto nel tuo ambiente
plot_data <- read.csv("D:/AcMet/3AcMet/Epigenetic_Analysis/Methylation/results/plot_data.csv", stringsAsFactors = FALSE)

cat("Dati letti da plot_data.csv\n")
cat("Dimensioni dataset:", dim(plot_data), "\n")
cat("Colonne:", names(plot_data), "\n")

# Controlla i valori unici di Gene e Project
cat("Geni presenti:", unique(plot_data$Gene), "\n")
cat("Progetti presenti:", unique(plot_data$Project), "\n")

# =========================================================================
# MODIFICA: Definisci l'ordine desiderato delle colonne dei progetti (tumori primari)
# Le colonne finali saranno: Probe, TCGA-LUAD, TCGA-PRAD, TCGA-COAD
project_col_order <- c("TCGA-LUAD", "TCGA-PRAD", "TCGA-COAD")
cat("\nOrdine colonne progetti desiderato:", project_col_order, "\n")
# =========================================================================

# Crea tabelle separate per CGAS e STING in formato wide
create_gene_tables <- function(plot_data, order_vec) {
  # Filtra e trasforma i dati per CGAS
  cgas_table <- plot_data %>%
    filter(Gene == "CGAS") %>%
    select(Probe, Project, BetaValue) %>%
    pivot_wider(
      names_from = Project,
      values_from = BetaValue,
      values_fn = mean  # Usa la media se ci sono duplicati
    ) %>%
    # =========================================================================
  # MODIFICA: Seleziona e riordina le colonne: Probe, seguite dall'ordine definito
  select(Probe, all_of(order_vec))
  # =========================================================================
  
  # Filtra e trasforma i dati per STING
  sting_table <- plot_data %>%
    filter(Gene == "STING") %>%
    select(Probe, Project, BetaValue) %>%
    pivot_wider(
      names_from = Project,
      values_from = BetaValue,
      values_fn = mean  # Usa la media se ci sono duplicati
    ) %>%
    # =========================================================================
  # MODIFICA: Seleziona e riordina le colonne: Probe, seguite dall'ordine definito
  select(Probe, all_of(order_vec))
  # =========================================================================
  
  return(list(CGAS = cgas_table, STING = sting_table))
}

# Crea le tabelle, passando il nuovo ordine di colonne
# =========================================================================
tables <- create_gene_tables(plot_data, project_col_order)
# =========================================================================

# Stampa anteprima
cat("\nAnteprima tabella CGAS (Verifica ordine colonne):\n")
print(head(tables$CGAS))
cat("\nAnteprima tabella STING (Verifica ordine colonne):\n")
print(head(tables$STING))

# Crea directory results se non esiste
if (!dir.exists("results")) {
  dir.create("results")
}

# Salva come file Excel
output_file <- "D:/AcMet/3AcMet/Epigenetic_Analysis/Methylation/results/methylation_beta_values_ordered.xlsx"
write_xlsx(tables, output_file)

cat("\n✓ File Excel creato con successo:", output_file, "\n")
cat("✓ Fogli creati: CGAS e STING\n")
cat("✓ Struttura: Righe = Sonde, Colonne = Tumori primari (LUAD, PRAD, COAD) - Ordine forzato\n")

# Verifica il file creato
cat("\nVerifica file:\n")
file_info <- file.info(output_file)
cat("Dimensione file:", round(file_info$size / 1024, 2), "KB\n")
cat("Path assoluto:", normalizePath(output_file), "\n")