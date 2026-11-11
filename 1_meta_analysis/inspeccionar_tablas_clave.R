# =====================================================================
# Script de Inspección Rápida para Tablas Clave del Pipeline 16S
# Ubicación: 1_meta_analysis/inspeccionar_tablas_clave.R
# USO: Rscript 1_meta_analysis/inspeccionar_tablas_clave.R
# =====================================================================

# Cargar setup global asegurando la ruta correcta desde la carpeta 1_meta_analysis/
if (file.exists("1_meta_analysis/00_setup.R")) {
  source("1_meta_analysis/00_setup.R")
} else if (file.exists("00_setup.R")) {
  source("00_setup.R")
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# Lista de las tablas principales (rutas relativas a la raíz del proyecto 16s_insomnia/)
archivos_clave <- c(
  "Diversidad Alfa"                  = "outputs/tables/alpha_diversity.csv",
  "Ratio Firmicutes/Bacteroidota"    = "outputs/tables/fb_ratio.csv",
  "ASVs Seleccionados (Boruta)"      = "outputs/boruta/boruta_selected_otus_with_taxonomy.csv",
  "Biomarcadores LEfSe"             = "outputs/lefse/lefse_sig_with_otu.csv",
  "Correlaciones Taxón-KO (Edges)"   = "outputs/function/tables/taxon_ko_edges_filtered.csv",
  "Fold-Change de KOs"               = "outputs/tables/KO_FC_group.csv"
)

cat("=====================================================================\n")
cat("INSPECCIÓN AUTOMATIZADA DE RESULTADOS CLAVE\n")
cat("=====================================================================\n\n")

for (nombre in names(archivos_clave)) {
  ruta <- archivos_clave[nombre]
  
  cat(sprintf("---------------------------------------------------------------------\n"))
  cat(sprintf("📁 TABLA: %s\n", nombre))
  cat(sprintf("📂 Archivo: %s\n", ruta))
  cat(sprintf("---------------------------------------------------------------------\n"))
  
  if (file.exists(ruta)) {
    # Cargar la tabla de forma segura
    df <- read_csv(ruta, show_col_types = FALSE)
    
    # Imprimir dimensiones y estructura básica
    cat(sprintf("Dimensiones: %d filas × %d columnas\n\n", nrow(df), ncol(df)))
    
    cat("Primeras 3 filas:\n")
    print(head(df, 3))
    cat("\n")
    
    # Mostrar un summary rápido
    if (nrow(df) > 0) {
      cat("Resumen estadístico básico:\n")
      print(summary(df))
    }
    
  } else {
    cat("⚠️ [AVISO] El archivo no se encontró en esta ruta (verifica si el módulo previo ya corrió).\n")
  }
  
  cat("\n\n")
}

cat("=====================================================================\n")
cat("¡Inspección finalizada con éxito!\n")
cat("=====================================================================\n")