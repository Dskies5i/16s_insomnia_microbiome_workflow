# ============================================================
# NOMBRE: 15_function_chord_from_ko_fc.R
# DESCRIPCIÓN: Gráfico de acorde basado en cambios de pliegue (FC) de KOs
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# USAGE: Rscript 1_meta_analysis/15_function_chord_from_ko_fc.R
# ============================================================

source("1_meta_analysis/00_setup.R")
source("utils/utils_io.R")
source("utils/utils_helper_functions.R")
library(circlize)
library(readr)
library(dplyr)
library(fs)

cat("Running Script 15: Functional Chord Diagram from KO Fold Change...\n")

out_dir <- file.path(OUT_DIR, "function")
dir_create(file.path(out_dir, "figures"))

fc_path <- "outputs/tables/KO_FC_group.csv"
if (!file.exists(fc_path)) {
  stop("X Error: KO_FC_group.csv no encontrado. Ejecuta primero el script 13a.")
}

fc_df <- read_csv(fc_path, show_col_types = FALSE)

if (nrow(fc_df) > 40) {
  message(sprintf("Limitando el diagrama de acorde a los primeros 40 KOs (de %d) para claridad visual.", nrow(fc_df)))
  fc_df <- fc_df %>% arrange(desc(abs(fc))) %>% head(40)
}

# FIX: Replaced syntax error %s% with standard %>% pipe
chord_df <- fc_df %>%
  mutate(
    value = abs(fc),
    direction = ifelse(fc > 0, "Insomnia", "Control"),
    edge_col = ifelse(fc > 0, "#d95f02", "#2b5c8f")
  ) %>%
  select(from = KO, to = direction, value, edge_col)

sectors <- unique(c(chord_df$from, chord_df$to))
grid_cols <- setNames(rep("#bdbdbd", length(sectors)), sectors)
grid_cols["Insomnia"] <- "#d95f02"
grid_cols["Control"]  <- "#2b5c8f"

png_path <- file.path(out_dir, "figures/chord_ko_fc.png")

png(png_path, width = 3600, height = 3600, res = 300)
circos.clear()
circos.par(start.degree = 90, gap.after = 5)

chordDiagram(
  x = chord_df[, c("from", "to", "value")],
  grid.col = grid_cols,
  col = chord_df$edge_col,
  transparency = 0.4,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    nm <- get.cell.meta.data("sector.index")
    xc <- get.cell.meta.data("xcenter")
    circos.text(xc, y = 0.5, labels = nm, facing = "clockwise", niceFacing = TRUE, cex = 0.6, font = 2)
  },
  bg.border = NA
)

circos.clear()
dev.off()

cat(sprintf("Gráfico de acorde generado con éxito en: %s\n", png_path))


























# 
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_io.R")
# source("utils/utils_helper_functions.R")
# 
# library(circlize)
# library(readr)
# library(dplyr)
# library(fs)
# 
# cat("Running Script 15: Functional Chord Diagram from KO Fold Change...\n")
# 
# # 1. Definir directorios de salida
# out_dir <- file.path(OUT_DIR, "function")
# dir_create(file.path(out_dir, "figures"))
# 
# # 2. Cargar datos de contribución diferencial o fold change de KOs
# fc_path <- "outputs/tables/KO_FC_group.csv"
# if (!file.exists(fc_path)) {
#   stop("❌ Error: KO_FC_group.csv no encontrado. Ejecuta primero el script 13a.")
# }
# 
# fc_df <- read_csv(fc_path, show_col_types = FALSE)
# 
# # 3. Preparar datos para el diagrama de acorde
# # Filtrar o seleccionar las características más relevantes si es necesario
# if (nrow(fc_df) > 40) {
#   message(sprintf("⚠️ Limitando el diagrama de acorde a los primeros 40 KOs (de %d) para claridad visual.", nrow(fc_df)))
#   fc_df <- fc_df %>% arrange(desc(abs(fc))) %>% head(40)
# }
# 
# # Construir tabla de enlaces/sectores para circlize
# chord_df <- fc_df %>%
#   mutate(
#     value = abs(fc),
#     direction = ifelse(fc > 0, "Insomnia", "Control"),
#     edge_col = ifelse(fc > 0, "#d95f02", "#2b5c8f")
#   ) %>%
#   select(from = KO, to = direction, value, edge_col)
# 
# # 4. Configurar paletas y colores del circlize
# sectors <- unique(c(chord_df$from, chord_df$to))
# grid_cols <- setNames(rep("#bdbdbd", length(sectors)), sectors)
# grid_cols["Insomnia"] <- "#d95f02"
# grid_cols["Control"] <- "#2b5c8f"
# 
# # 5. Generar PDF y PNG del diagrama de acorde
# pdf_path <- file.path(out_dir, "figures/chord_ko_fc.pdf")
# png_path <- file.path(out_dir, "figures/chord_ko_fc.png")
# 
# # PNG de alta resolución
# png(png_path, width = 3600, height = 3600, res = 300)
# circos.clear()
# circos.par(start.degree = 90, gap.after = 5)
# chordDiagram(
#   x = chord_df[, c("from", "to", "value")],
#   grid.col = grid_cols,
#   col = chord_df$edge_col,
#   transparency = 0.4,
#   annotationTrack = "grid",
#   preAllocateTracks = list(track.height = 0.1)
# )
# circos.trackPlotRegion(
#   track.index = 1, 
#   panel.fun = function(x, y) {
#     nm <- get.cell.meta.data("sector.index")
#     xc <- get.cell.meta.data("xcenter")
#     circos.text(xc, y = 0.5, labels = nm, facing = "clockwise", niceFacing = TRUE, cex = 0.6, font = 2)
#   }, 
#   bg.border = NA
# )
# circos.clear()
# dev.off()
# 
# cat(sprintf("✅ Gráfico de acorde generado con éxito en: %s\n", png_path))
# 