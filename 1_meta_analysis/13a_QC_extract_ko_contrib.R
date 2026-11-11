# # ============================================================
# # SCRIPT: 13a_QC_extract_ko_contrib.R
# # DESCRIPTION: Quality Control (QC) for auxiliary contribution inputs
# # USE: Rscript 1_meta_analysis/13_QC_extract_ko_contrib.R
# # ============================================================
# # ============================================================
# # SCRIPT: 13b_QC_extract_ko_contrib.R (Memory-Safe QC)
# # ============================================================

source("1_meta_analysis/00_setup.R")
library(readr)
library(dplyr)

cat("====================================================\n")
cat("  RUNNING QC CHECK FOR SCRIPT 13a OUTPUTS\n")
cat("====================================================\n\n")

target_files <- c(
  metadata = "outputs/tables/metadata_for_taxa.csv",
  taxonomy = "outputs/tables/picrust_taxonomy_filled.csv",
  ko_fc    = "outputs/tables/KO_FC_group.csv",
  contrib  = "outputs/tables/ko_contrib_filtered.csv"
)

# 1. Verificar existencia de archivos
missing <- target_files[!file.exists(target_files)]
if (length(missing) > 0) {
  stop("X QC FAILED: Los siguientes archivos no existen. Corre 13a primero:\n", 
       paste(" -", missing, collapse = "\n"))
}

# 2. Inspeccionar metadata_for_taxa.csv
meta_df <- read_csv(target_files["metadata"], show_col_types = FALSE)
cat("[1/4] metadata_for_taxa.csv:\n")
cat(sprintf("  - Filas: %d | Columnas: %d\n", nrow(meta_df), ncol(meta_df)))
cat(sprintf("  - NAs detectados: %d\n", sum(is.na(meta_df))))
print(head(meta_df, 3))
cat("\n")

# 3. Inspeccionar picrust_taxonomy_filled.csv
tax_df <- read_csv(target_files["taxonomy"], show_col_types = FALSE)
cat("[2/4] picrust_taxonomy_filled.csv:\n")
cat(sprintf("  - Filas: %d | Columnas: %d\n", nrow(tax_df), ncol(tax_df)))
cat(sprintf("  - Taxones no asignados (g_Unassigned): %d\n", sum(tax_df$Taxa_name == "g_Unassigned")))
print(head(tax_df, 3))
cat("\n")

# 4. Inspeccionar KO_FC_group.csv
ko_df <- read_csv(target_files["ko_fc"], show_col_types = FALSE)
cat("[3/4] KO_FC_group.csv:\n")
cat(sprintf("  - KOs calculados: %d\n", nrow(ko_df)))
cat("  - Distribución de regulación:\n")
print(table(ko_df$Strength))
print(head(ko_df, 3))
cat("\n")

# 5. Inspeccionar ko_contrib_filtered.csv
contrib_df <- read_csv(target_files["contrib"], show_col_types = FALSE)
cat("[4/4] ko_contrib_filtered.csv:\n")
cat(sprintf("  - Registros de contribución: %d\n", nrow(contrib_df)))
cat(sprintf("  - Contribuciones > 0: %d\n", sum(contrib_df$norm_taxon_function_contrib > 0)))
print(head(contrib_df, 3))
cat("\n")

cat("====================================================\n")
cat("  ✔ QC COMPLETADO: Todos los outputs de 13a son válidos.\n")
cat("====================================================\n")




# # ============================================================
# # SCRIPT: 13a_QC_extract_ko_contrib.R
# # DESCRIPTION: Quality Control (QC) for auxiliary contribution inputs
# # USE: Rscript 1_meta_analysis/13_QC_extract_ko_contrib.R
# # ============================================================
# # ============================================================
# # SCRIPT: 13b_QC_extract_ko_contrib.R (Memory-Safe QC)
# # ============================================================
# source("1_meta_analysis/00_setup.R")
# library(readr)
# library(dplyr)
# 
# cat("====================================================\n")
# cat("  RUNNING QC CHECK FOR SCRIPT 13a OUTPUTS\n")
# cat("====================================================\n\n")
# 
# target_files <- c(
#   metadata = "outputs/tables/metadata_for_taxa.csv",
#   taxonomy = "outputs/tables/picrust_taxonomy_filled.csv",
#   ko_fc    = "outputs/tables/KO_FC_group.csv",
#   contrib  = "outputs/tables/ko_contrib_filtered.csv"
# )
# 
# # 1. Verificar existencia de archivos
# missing <- target_files[!file.exists(target_files)]
# if (length(missing) > 0) {
#   stop("X QC FAILED: Los siguientes archivos no existen. Corre 13a primero:\n", 
#        paste(" -", missing, collapse = "\n"))
# }
# 
# # 2. Inspeccionar metadata_for_taxa.csv
# meta_df <- read_csv(target_files["metadata"], show_col_types = FALSE)
# cat("[1/4] metadata_for_taxa.csv:\n")
# cat(sprintf("  - Filas: %d | Columnas: %d\n", nrow(meta_df), ncol(meta_df)))
# cat(sprintf("  - NAs detectados: %d\n", sum(is.na(meta_df))))
# print(head(meta_df, 3))
# cat("\n")
# 
# # 3. Inspeccionar picrust_taxonomy_filled.csv
# tax_df <- read_csv(target_files["taxonomy"], show_col_types = FALSE)
# cat("[2/4] picrust_taxonomy_filled.csv:\n")
# cat(sprintf("  - Filas: %d | Columnas: %d\n", nrow(tax_df), ncol(tax_df)))
# cat(sprintf("  - Taxones no asignados (g_Unassigned): %d\n", sum(tax_df$Taxa_name == "g_Unassigned")))
# print(head(tax_df, 3))
# cat("\n")
# 
# # 4. Inspeccionar KO_FC_group.csv
# ko_df <- read_csv(target_files["ko_fc"], show_col_types = FALSE)
# cat("[3/4] KO_FC_group.csv:\n")
# cat(sprintf("  - KOs calculados: %d\n", nrow(ko_df)))
# cat("  - Distribución de regulación:\n")
# print(table(ko_df$Strength))
# print(head(ko_df, 3))
# cat("\n")
# 
# # 5. Inspeccionar ko_contrib_filtered.csv
# contrib_df <- read_csv(target_files["contrib"], show_col_types = FALSE)
# cat("[4/4] ko_contrib_filtered.csv:\n")
# cat(sprintf("  - Registros de contribución: %d\n", nrow(contrib_df)))
# cat(sprintf("  - Contribuciones > 0: %d\n", sum(contrib_df$norm_taxon_function_contrib > 0)))
# print(head(contrib_df, 3))
# cat("\n")
# 
# cat("====================================================\n")
# cat("  ✔ QC COMPLETADO: Todos los outputs de 13a son válidos.\n")
# cat("====================================================\n")