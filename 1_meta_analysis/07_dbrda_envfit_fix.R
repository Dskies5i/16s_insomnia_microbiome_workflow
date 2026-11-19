##NEWER NEWER VERSION IN LINUX TOO
# ============================================================
# NOMBRE: 07_dbrda_envfit_fix.R
# DESCRIPCIÓN: Análisis db-RDA y selección de modelos basados en 
#              distancias con ajuste de covariables, envfit taxonómico
#              y generación de gráficos descriptivos (PDF + PNG).
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# USAGE: Rscript 1_meta_analysis/07_dbrda_envfit_fix.R
# ============================================================

source("1_meta_analysis/00_setup.R")
source("utils/utils_phyloseq.R")
source("utils/utils_io.R")
source("utils/utils_plot.R")

library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(ggplot2)

# 1. Primero cargamos la configuración para que 'cfg' esté disponible
cfg <- yaml::read_yaml("config/config.yaml")

out_dir <- file.path(cfg$outputs$dir, "dbrda")
dir_create(out_dir)

# ============================================================
# CHECKPOINT: Omitir computación pesada si ya existen los outputs
# ============================================================
final_target_file <- file.path(out_dir, "dbrda_model.rds")

if (file.exists(final_target_file)) {
  cat(sprintf("⏩ [SALTADO] Los resultados de db-RDA ya existen en: %s\n", final_target_file))
  cat("   (Borra 'outputs/dbrda/dbrda_model.rds' si deseas regenerar modelos y gráficos).\n\n")
} else {
  
  # ============================================================
  # SELECCIÓN DE MODELO: FILTRADO VS UNFILTERED
  # (Descomenta la opción que prefieras usar)
  # ============================================================
  
  # OPCIÓN A: Usar el objeto filtrado por defecto (Alineado con Boruta/LEfSe/MaAsLin2)
  physeq <- readRDS(file.path(cfg$outputs$dir, "models/physeq_filtered.rds"))
  model_type_used <- "FILTERED"
  
  # OPCIÓN B: Usar el objeto crudo/unfiltered (Comenta la línea de arriba y descomenta esta)
  # physeq <- build_phyloseq()
  # model_type_used <- "UNFILTERED"
  
  pseq <- physeq
  cat(sprintf("✅ RUNNING %s MODEL.\n", model_type_used))
  
  dist_method <- cfg$analysis$beta$distance
  bray_dist <- phyloseq::distance(pseq, method = dist_method)
  
  # Convertir explícitamente a data.frame plano para compatibilidad con vegan
  meta <- as.data.frame(as(sample_data(pseq), "data.frame"), stringsAsFactors = FALSE)
  group_col <- cfg$columns$group
  meta[[group_col]] <- factor(meta[[group_col]], levels = cfg$groups$levels)
  
  # Sincronizar muestras con labels en minúscula por seguridad
  common_samples <- intersect(rownames(meta), labels(bray_dist))
  meta <- meta[common_samples, , drop = FALSE]
  bray_dist <- as.dist(as.matrix(bray_dist)[common_samples, common_samples])
  
  # ============================================================
  # SELECCIÓN ESCALONADA (STEPWISE SCOPE)
  # ============================================================
  scope_covs <- cfg$differential$dbrda$scope_covariates
  mod0 <- stats::as.formula(paste("bray_dist ~", group_col))
  mod1 <- stats::as.formula(paste("bray_dist ~", group_col, "+", paste(scope_covs, collapse = " + ")))
  
  cap0 <- capscale(mod0, data = meta)
  cap1 <- capscale(mod1, data = meta)
  
  step_model <- ordiR2step(cap0, scope = formula(cap1),
                           direction = "forward", R2scope = TRUE, 
                           pstep = cfg$differential$dbrda$permutations)
  
  saveRDS(step_model, file.path(out_dir, "ordiR2step_model.rds"))
  writeLines(capture.output(step_model), con = file.path(out_dir, "ordiR2step_log.txt"))
  
  # ============================================================
  # MODELO CONDICIONADO db-RDA Y ESTADÍSTICAS
  # ============================================================
  cond_covs <- cfg$differential$dbrda$condition_covariates
  if (length(cond_covs) > 0) {
    dbRDA_formula <- stats::as.formula(paste0("bray_dist ~ ", group_col, " + Condition(", paste(cond_covs, collapse = ", "), ")"))
    cond_text <- paste(cond_covs, collapse = ", ")
  } else {
    dbRDA_formula <- stats::as.formula(paste0("bray_dist ~ ", group_col))
    cond_text <- "None"
  }
  
  dbRDA_model <- capscale(dbRDA_formula, data = meta)
  saveRDS(dbRDA_model, file.path(out_dir, "dbrda_model.rds"))
  
  anova_all <- anova(dbRDA_model, permutations = cfg$differential$dbrda$permutations)
  anova_margin <- anova(dbRDA_model, by = "margin", permutations = cfg$differential$dbrda$permutations)
  
  write_csv_safe(as.data.frame(anova_all) %>% tibble::rownames_to_column("Term"), file.path(out_dir, "dbrda_anova_all.csv"))
  write_csv_safe(as.data.frame(anova_margin) %>% tibble::rownames_to_column("Term"), file.path(out_dir, "dbrda_anova_margin.csv"))
  
  # R² y R² Ajustado
  r2 <- RsquareAdj(dbRDA_model)
  writeLines(c(
    paste0("R2:           ", r2$r.squared),
    paste0("Adjusted R2:  ", r2$adj.r.squared)
  ), con = file.path(out_dir, "dbrda_r2.txt"))
  
  # ============================================================
  # ENVFIT: FILTRO DE PREVALENCIA Y VECTORES TAXONÓMICOS
  # ============================================================
  otu_mat <- as(otu_table(pseq), "matrix")
  if (!taxa_are_rows(pseq)) otu_mat <- t(otu_mat)
  otu_mat <- t(otu_mat)[common_samples, , drop = FALSE]
  
  prev <- apply(otu_mat, 2, function(x) sum(x > 0) / length(x))
  otu_f <- otu_mat[, prev > cfg$differential$dbrda$prevalence_min, drop = FALSE]
  
  fit <- envfit(dbRDA_model, otu_f,
                perm  = cfg$differential$dbrda$permutations,
                p.max = 1)
  
  fit_scores <- as.data.frame(scores(fit, display = "vectors"))
  fit_scores$pval <- fit$vectors$pvals
  fit_scores$r2_vec <- fit$vectors$r
  fit_scores$OTU_ID <- rownames(fit_scores)
  
  stopifnot(nrow(fit_scores) == length(fit$vectors$pvals))
  
  fit_sig <- fit_scores[!is.na(fit_scores$pval) & fit_scores$pval < 0.05, , drop = FALSE]
  
  tax_df <- as.data.frame(tax_table(pseq)) %>% tibble::rownames_to_column("OTU_ID")
  fit_sig_annot <- merge(fit_sig, tax_df, by = "OTU_ID") %>%
    arrange(pval, desc(r2_vec))
  
  write_csv_safe(fit_sig, file.path(out_dir, "envfit_significant_otus.csv"))
  write_csv_safe(fit_sig_annot, file.path(out_dir, "envfit_significant_otus_with_taxonomy.csv"))
  
  # ============================================================
  # EXTRACCIÓN DE MÉTRICAS PARA VISUALIZACIÓN DESCRIPTIVA
  # ============================================================
  sites_df <- as.data.frame(scores(dbRDA_model, display = "sites", choices = c(1, 2)))
  axis_names <- colnames(sites_df) 
  colnames(sites_df) <- c("Axis1", "Axis2")
  sites_df$Sample <- rownames(sites_df)
  sites_df$Group <- meta[sites_df$Sample, group_col]
  
  eig_summary <- summary(dbRDA_model)$cont$importance
  var_ax1 <- round(eig_summary["Proportion Explained", axis_names[1]] * 100, 2)
  var_ax2 <- round(eig_summary["Proportion Explained", axis_names[2]] * 100, 2)
  
  anova_df <- as.data.frame(anova_margin)
  f_val <- round(anova_df[1, "F"], 3)
  p_val <- anova_df[1, "Pr(>F)"]
  p_text <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.4f", p_val))
  
  r2_pct <- round(r2$r.squared * 100, 2)
  adj_r2_pct <- round(r2$adj.r.squared * 100, 2)
  
  group_counts <- table(sites_df$Group)
  counts_str <- paste(names(group_counts), group_counts, sep = ": n=", collapse = " | ")
  
  stats_box <- sprintf(
    "Modelo: ~ %s + Condition(%s)\nDistancia: %s | Permutaciones: %d\nR² = %.2f%% (Adj. R² = %.2f%%)\nANOVA (%s): F = %.3f, %s\nASVs sig. envfit (p < 0.05): %d / %d",
    group_col, cond_text, toupper(dist_method), cfg$differential$dbrda$permutations,
    r2_pct, adj_r2_pct, group_col, f_val, p_text, nrow(fit_sig_annot), ncol(otu_f)
  )
  
  colors_map <- unlist(cfg$groups$colors)
  
  # ============================================================
  # GRÁFICO 1: ORDENACIÓN DE MUESTRAS (SITES) CON ESTADÍSTICAS
  # ============================================================
  p_sites <- ggplot(sites_df, aes(x = Axis1, y = Axis2, color = Group, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, linewidth = 0.6, show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(values = colors_map) +
    scale_fill_manual(values = colors_map) +
    annotate(
      "label",
      x = -Inf, y = Inf,
      label = stats_box,
      hjust = -0.04, vjust = 1.08,
      size = 3.1,
      fill = alpha("white", 0.88),
      color = "grey15",
      label.size = 0.3
    ) +
    labs(
      title = "Distance-based Redundancy Analysis (db-RDA)",
      subtitle = sprintf("Muestras totales: N = %d (%s) | Elipses al 95%% de confianza", nrow(sites_df), counts_str),
      x = sprintf("%s (%.2f%% varianza explicada / constrained)", axis_names[1], var_ax1),
      y = sprintf("%s (%.2f%% varianza residual / unconstrained)", axis_names[2], var_ax2),
      color = "Grupo",
      caption = "Generado por Script 07 (vegan::capscale) | Covariables ajustadas mediante Condition()"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(out_dir, "dbrda_sites.pdf"), plot = p_sites, width = 8.5, height = 6.5)
  ggsave(file.path(out_dir, "dbrda_sites.png"), plot = p_sites, width = 8.5, height = 6.5, dpi = 300)
  
  # ============================================================
  # GRÁFICO 2 (BONUS): BIPLOT CON TOP 10 VECTORES ENVFIT
  # ============================================================
  if (nrow(fit_sig_annot) > 0) {
    top_vec <- head(fit_sig_annot, 10)
    colnames(top_vec)[colnames(top_vec) == axis_names[1]] <- "Vec1"
    colnames(top_vec)[colnames(top_vec) == axis_names[2]] <- "Vec2"
    
    top_vec$TaxLabel <- ifelse(
      !is.na(top_vec$Genus) & top_vec$Genus != "" & top_vec$Genus != "Unassigned",
      paste0(gsub("^g__", "", top_vec$Genus), " (", substr(top_vec$OTU_ID, 1, 5), ")"),
      paste0("ASV_", substr(top_vec$OTU_ID, 1, 6))
    )
    
    mult <- min(
      (max(sites_df$Axis1) - min(sites_df$Axis1)) / (max(abs(top_vec$Vec1)) * 2.2),
      (max(sites_df$Axis2) - min(sites_df$Axis2)) / (max(abs(top_vec$Vec2)) * 2.2)
    )
    if (!is.finite(mult) || mult <= 0) mult <- 1
    
    top_vec$xend <- top_vec$Vec1 * mult
    top_vec$yend <- top_vec$Vec2 * mult
    
    p_biplot <- p_sites +
      geom_segment(
        data = top_vec,
        aes(x = 0, y = 0, xend = xend, yend = yend),
        inherit.aes = FALSE,
        arrow = arrow(length = unit(0.22, "cm")),
        color = "#1b7837",
        linewidth = 0.65,
        alpha = 0.85
      ) +
      geom_text(
        data = top_vec,
        aes(x = xend * 1.08, y = yend * 1.08, label = TaxLabel),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "italic",
        color = "#00441b"
      ) +
      labs(
        title = "db-RDA Biplot con Vectores Taxonómicos Significativos (envfit)",
        subtitle = sprintf("Top %d ASVs asociados (p < 0.05) | R² = %.2f%% (%s)", nrow(top_vec), r2_pct, p_text)
      )
    
    ggsave(file.path(out_dir, "dbrda_biplot_envfit.pdf"), plot = p_biplot, width = 9, height = 7)
    ggsave(file.path(out_dir, "dbrda_biplot_envfit.png"), plot = p_biplot, width = 9, height = 7, dpi = 300)
  }
  
  cat("✅ Script 07 finalizado con éxito. Gráficos guardados en PDF y PNG en outputs/dbrda/.\n")
}



# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================



# 
# # ============================================================
# # NOMBRE: 07_dbrda_envfit_fix.R
# # DESCRIPCIÓN: Análisis db-RDA y selección de modelos basados en 
# #              distancias con ajuste de covariables y envfit taxonómico.
# # 
# # ¿QUÉ HACE ESTE SCRIPT?
# #    1. Carga los datos usando las funciones estándar del pipeline.
# #    2. Ejecuta selección de modelos escalonados (ordiR2step).
# #    3. Ajusta el modelo db-RDA final por condición y calcula pseudo-R2 y ANOVA.
# #    4. Aplica el análisis envfit original sobre los ASVs filtrados.
# #
# # ARCHIVOS DE ENTRADA (INPUTS):
# #    - outputs/models/physeq_filtered.rds                    (Objeto phyloseq filtrado)
# #    - config/config.yaml                                    (Archivo de configuración)
# #
# # USO: Rscript 1_meta_analysis/07_dbrda_envfit_fix.R
# # ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# # ============================================================
# 
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_phyloseq.R")
# source("utils/utils_io.R")
# 
# # 1. Primero cargamos la configuración para que 'cfg' esté disponible
# cfg <- yaml::read_yaml("config/config.yaml")
# 
# out_dir <- file.path(cfg$outputs$dir, "dbrda")
# dir_create(out_dir)
# 
# # ============================================================
# # CHECKPOINT: Omitir computación pesada si ya existen los outputs
# # ============================================================
# final_target_file <- file.path(out_dir, "dbrda_model.rds")
# 
# if (file.exists(final_target_file)) {
#   cat(sprintf("⏩ [SALTADO] Los resultados de db-RDA ya existen en: %s\n", final_target_file))
#   cat("   (Se omite la computación pesada y se continúa con el pipeline).\n\n")
# } else {
#   
#   # ============================================================
#   # CARGA DE DATOS Y PREPARACIÓN
#   #
#   ###=========FILTERED VS UNFILTERED=======================
#   # ============================================================
#   ###==========================================================
#   ###=========================================================
# 
#   ## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#   
#   ## OPTION A: FILTERED (DEFAULT) ================
#   # CHANGE A: Usamos el objeto phyloseq FILTRADO (2394 taxa tras
#   # g__other|Unassigned + prune cero).  Esto alinea 07 con el resto
#   # del pipeline (02–06, 08, 09) — misma feature universe, permite
#   # cruzar hits de Boruta/LEfSe con envfit.  La prevalencia (>5%)
#   # se aplica igualmente abajo como filtro estadístico adicional.
#   # The filtered model (outputs/models/physeq_filtered.rds) is generated by the script 1_meta_analysis/01_build_phyloseq.R
#   
#   physeq <- readRDS(file.path(cfg$outputs$dir, "models/physeq_filtered.rds"))
#   pseq <- physeq
#   cat("✅ RUNNING FILTERED MODEL.\n")
# 
#   ## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#   
#   ### OPTION B UNFILTERED (DESACTIVADO) ==========
#   # NOTE: We deliberately use the RAW phyloseq (not physeq_filtered.rds)
#   # because db-RDA+envfit benefits from including g__other taxa for
#   # full community analysis. Prevalence filtering (5%) is applied below
#   # as a statistical filter, not taxonomic. This differs from Boruta
#   # and LEfSe which use physeq_filtered.rds (g__other-filtered).
#   
#   #physeq <- build_phyloseq()
#   #pseq <- physeq
#   #cat("✅ RUNNING UNFILTERED MODEL.\n")
# 
#   ## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#   
#   # ============================================================
#   # NOTA TÉCNICA: MODELO FILTRADO VS. UNFILTERED EN DB-RDA / ENVFIT
#   # ============================================================
#   # 
#   # 1. MODELO FILTRADO (DEFAULT - physeq_filtered.rds):
#   #    - Elimina taxones no deseados como 'g__other' y 'Unassigned'[cite: 15].
#   #    - Mantiene la consistencia con el resto del pipeline analítico 
#   #      (02–06, 08, 09) usando exactamente el mismo "feature universe"[cite: 10].
#   #    - Facilita el cruce y validación cruzada entre los hits de envfit 
#   #      y los biomarcadores de Boruta o LEfSe[cite: 10].
#   #
#   # 2. MODELO UNFILTERED (DESACTIVADO - physeq crudo / build_phyloseq()):
#   #    - Conserva la totalidad de la comunidad microbiana, incluyendo 
#   #      categorías no asignadas y taxones generales 'g__other'[cite: 10].
#   #    - Útil para análisis comunitarios globales de ordenación (db-RDA), 
#   #      aplicando el filtro de prevalencia (>5%) de forma puramente estadística[cite: 10].
#   #    - Crea una divergencia de features frente a Boruta y LEfSe (que usan 
#   #      el objeto filtrado), lo que dificulta comparar directamente los resultados[cite: 10].
#   # ============================================================
# 
#   
#   dist_method <- cfg$analysis$beta$distance
#   bray_dist <- phyloseq::distance(pseq, method = dist_method)
#   
#   # Convertir explícitamente a data.frame plano para compatibilidad con vegan
#   meta <- as.data.frame(as(sample_data(pseq), "data.frame"), stringsAsFactors = FALSE)
#   group_col <- cfg$columns$group
#   meta[[group_col]] <- factor(meta[[group_col]], levels = cfg$groups$levels)
#   
#   # Sincronizar muestras con labels en minúscula por seguridad
#   common_samples <- intersect(rownames(meta), labels(bray_dist))
#   meta <- meta[common_samples, , drop = FALSE]
#   bray_dist <- as.dist(as.matrix(bray_dist)[common_samples, common_samples])
#   
#   # stepwise scope
#   scope_covs <- cfg$differential$dbrda$scope_covariates
#   mod0 <- stats::as.formula(paste("bray_dist ~", group_col))
#   mod1 <- stats::as.formula(paste("bray_dist ~", group_col, "+", paste(scope_covs, collapse = " + ")))
#   
#   cap0 <- capscale(mod0, data = meta)
#   cap1 <- capscale(mod1, data = meta)
#   
#   step_model <- ordiR2step(cap0, scope = formula(cap1),
#                            direction = "forward", R2scope = TRUE, pstep = cfg$differential$dbrda$permutations)
#   
#   saveRDS(step_model, file.path(out_dir, "ordiR2step_model.rds"))
#   writeLines(capture.output(step_model), con = file.path(out_dir, "ordiR2step_log.txt"))
#   
#   # Condition modelo
#   cond_covs <- cfg$differential$dbrda$condition_covariates
#   if (length(cond_covs) > 0) {
#     dbRDA_formula <- stats::as.formula(paste0("bray_dist ~ ", group_col, " + Condition(", paste(cond_covs, collapse = ", "), ")"))
#   } else {
#     dbRDA_formula <- stats::as.formula(paste0("bray_dist ~ ", group_col))
#   }
#   
#   dbRDA_model <- capscale(dbRDA_formula, data = meta)
#   saveRDS(dbRDA_model, file.path(out_dir, "dbrda_model.rds"))
#   
#   anova_all <- anova(dbRDA_model, permutations = cfg$differential$dbrda$permutations)
#   anova_margin <- anova(dbRDA_model, by = "margin", permutations = cfg$differential$dbrda$permutations)
#   
#   write_csv_safe(as.data.frame(anova_all) %>% tibble::rownames_to_column("Term"), file.path(out_dir, "dbrda_anova_all.csv"))
#   write_csv_safe(as.data.frame(anova_margin) %>% tibble::rownames_to_column("Term"), file.path(out_dir, "dbrda_anova_margin.csv"))
#   
#   # R²
#   r2 <- RsquareAdj(dbRDA_model)
#   writeLines(c(
#     paste0("R2:           ", r2$r.squared),
#     paste0("Adjusted R2:  ", r2$adj.r.squared)
#   ), con = file.path(out_dir, "dbrda_r2.txt"))
#   
#   # envfit: prevalence filter
#   otu_mat <- as(otu_table(pseq), "matrix")
#   if (!taxa_are_rows(pseq)) otu_mat <- t(otu_mat)
#   otu_mat <- t(otu_mat)[common_samples, , drop = FALSE]
#   
#   prev <- apply(otu_mat, 2, function(x) sum(x > 0) / length(x))
#   otu_f <- otu_mat[, prev > cfg$differential$dbrda$prevalence_min, drop = FALSE]
#   
#   fit <- envfit(dbRDA_model, otu_f,
#                 perm  = cfg$differential$dbrda$permutations,
#                 p.max = 1)
#   
#   fit_scores <- as.data.frame(scores(fit, display = "vectors"))
#   fit_scores$pval <- fit$vectors$pvals
#   fit_scores$OTU_ID <- rownames(fit_scores)
#   
#   stopifnot(nrow(fit_scores) == length(fit$vectors$pvals))
#   
#   fit_sig <- fit_scores[!is.na(fit_scores$pval) & fit_scores$pval < 0.05, , drop = FALSE]
#   
#   tax_df <- as.data.frame(tax_table(pseq)) %>% tibble::rownames_to_column("OTU_ID")
#   fit_sig_annot <- merge(fit_sig, tax_df, by = "OTU_ID")
#   
#   write_csv_safe(fit_sig, file.path(out_dir, "envfit_significant_otus.csv"))
#   write_csv_safe(fit_sig_annot, file.path(out_dir, "envfit_significant_otus_with_taxonomy.csv"))
#   
#   # Gráfico de ordenación
#   pdf(file.path(out_dir, "dbrda_sites.pdf"), width = 6, height = 5)
#   plot(dbRDA_model, display = "sites", main = "db-RDA")
#   dev.off()
#   
#   cat("✅ Script 07 finalizado con éxito.\n")
# }




# ============================================================
# ============================================================
# ============================================================
# ============================================================
# ============================================================





## PREVIOUS WORKING VERSION WITHOUT SKIP-CHECK (KEEP FOR SECURITY):
# # ============================================================
# # NOMBRE: 07_dbrda_envfit_fix.R
# # DESCRIPCIÓN: Análisis db-RDA y selección de modelos basados en 
# #              distancias con ajuste de covariables y envfit taxonómico.
# # 
# # ¿QUÉ HACE ESTE SCRIPT?
# #    1. Carga los datos usando las funciones estándar del pipeline.
# #    2. Ejecuta selección de modelos escalonados (ordiR2step).
# #    3. Ajusta el modelo db-RDA final por condición y calcula pseudo-R2 y ANOVA.
# #    4. Aplica el análisis envfit original sobre los ASVs filtrados.
# #
# # ARCHIVOS DE ENTRADA (INPUTS):
# #    - outputs/models/physeq_filtered.rds                  (Objeto phyloseq filtrado)
# #    - config/config.yaml                                  (Archivo de configuración)
# #
# # ============================================================
# # NOTAS Y SALIDAS ESPERADAS DE LA EJECUCIÓN:
# # ============================================================
# #
# # 1. MENSAJE "null device 1":
# #    Aparece al ejecutar `dev.off()` tras generar el PDF (`dbrda_sites.pdf`).
# #    Indica que el dispositivo gráfico se cerró y guardó correctamente 
# #    retornando al estado por defecto ("null device"). Es un comportamiento normal.
# #
# # 2. ARCHIVOS GENERADOS EN 'outputs/dbrda/':
# #    - ordiR2step_model.rds                  : Objeto R con el modelo de selección stepwise.
# #    - ordiR2step_log.txt                    : Resumen en texto del proceso de selección de variables.
# #    - dbrda_model.rds                       : Modelo db-RDA final ajustado.
# #    - dbrda_anova_all.csv                   : Resultados de ANOVA global del modelo.
# #    - dbrda_anova_margin.csv                : Resultados de ANOVA marginal por variable.
# #    - dbrda_r2.txt                          : Valores de R2 y R2 ajustado del modelo final.
# #    - envfit_significant_otus.csv           : ASVs con p < 0.05 según las permutaciones de envfit.
# #    - envfit_significant_otus_with_taxonomy.csv : ASVs significativos anotados con su taxonomía.
# #    - dbrda_sites.pdf                       : Gráfico de ordenación db-RDA generado y cerrado con dev.off().
# # ============================================================
# #
# #
# # USO: Rscript 1_meta_analysis/07_dbrda_envfit_fix.R
# #
# # TARDA MUCHO EN 8 GB DE RAM
# #
# # ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# # ============================================================
# ############ WORKING VERSION (TOOK HOURS!!)
# 
# source("1_meta_analysis/00_setup.R")
# source("utils/utils_phyloseq.R")
# source("utils/utils_io.R")
# 
# # 1. Primero cargamos la configuración para que 'cfg' esté disponible
# cfg <- yaml::read_yaml("config/config.yaml")
# 
# out_dir <- file.path(cfg$outputs$dir, "dbrda")
# dir_create(out_dir)
# 
# # ============================================================
# # CHECKPOINT: Omitir si ya existen los outputs principales
# # ============================================================
# final_target_file <- file.path(out_dir, "dbrda_model.rds")
# if (file.exists(final_target_file)) {
#   cat(sprintf("⏩ [SALTADO] Los resultados de db-RDA ya existen en: %s\n", final_target_file))
#   cat("   (Borra los archivos de 'outputs/dbrda/' manualmente si deseas volver a ejecutarlo).\n")
#   quit(save = "no", status = 0)
# }
# 
# ###==========================================================
# ###=========FILTERED VS UNFILTERED=======================
# ###==========================================================
# ##=============== FILTERED (DEFAULT) ================
# # CHANGE A: Usamos el objeto phyloseq FILTRADO (2394 taxa tras
# # g__other|Unassigned + prune cero).  Esto alinea 07 con el resto
# # del pipeline (02–06, 08, 09) — misma feature universe, permite
# # cruzar hits de Boruta/LEfSe con envfit.  La prevalencia (>5%)
# # se aplica igualmente abajo como filtro estadístico adicional.
# 
# 
# physeq <- readRDS(file.path(cfg$outputs$dir, "models/physeq_filtered.rds"))
# pseq <- physeq
# cat("✅ RUNNING FILTERED MODEL.\n")
# 
# ###============== UNFILTERED (DESACTIVADO) ==========
# # NOTE: We deliberately use the RAW phyloseq (not physeq_filtered.rds)
# # because db-RDA+envfit benefits from including g__other taxa for 
# # full community analysis. Prevalence filtering (5%) is applied below
# # as a statistical filter, not taxonomic. This differs from Boruta
# # and LEfSe which use physeq_filtered.rds (g__other-filtered).
# 
# #physeq <- build_phyloseq()
# #pseq <- physeq
# #cat("✅ RUNNING UNFILTERED MODEL.\n")
# 
#  
# ###==========================================================
# ###==========================================================
# # Change this line:
# # dist_method <- cfg$differential$dbrda$distance
# 
# # To this:
# dist_method <- cfg$analysis$beta$distance
# bray_dist <- phyloseq::distance(pseq, method = dist_method)
# 
# # Convertir explícitamente a data.frame plano para compatibilidad con vegan
# meta <- as.data.frame(as(sample_data(pseq), "data.frame"), stringsAsFactors = FALSE)
# group_col <- cfg$columns$group
# meta[[group_col]] <- factor(meta[[group_col]], levels = cfg$groups$levels)
# 
# # Sincronizar muestras con labels en minúscula por seguridad
# common_samples <- intersect(rownames(meta), labels(bray_dist))
# meta <- meta[common_samples, , drop = FALSE]
# bray_dist <- as.dist(as.matrix(bray_dist)[common_samples, common_samples])
# 
# # stepwise scope
# scope_covs <- cfg$differential$dbrda$scope_covariates
# mod0 <- stats::as.formula(paste("bray_dist ~", group_col))
# mod1 <- stats::as.formula(paste("bray_dist ~", group_col, "+", paste(scope_covs, collapse = " + ")))
# 
# cap0 <- capscale(mod0, data = meta)
# cap1 <- capscale(mod1, data = meta)
# 
# step_model <- ordiR2step(cap0, scope = formula(cap1),
#                          direction = "forward", R2scope = TRUE, pstep = cfg$differential$dbrda$permutations)
# 
# saveRDS(step_model, file.path(out_dir, "ordiR2step_model.rds"))
# writeLines(capture.output(step_model), con = file.path(out_dir, "ordiR2step_log.txt"))
# 
# # Condition modelo
# cond_covs <- cfg$differential$dbrda$condition_covariates
# if (length(cond_covs) > 0) {
#   dbRDA_formula <- stats::as.formula(paste0("bray_dist ~ ", group_col, " + Condition(", paste(cond_covs, collapse = ", "), ")"))
# } else {
#   dbRDA_formula <- stats::as.formula(paste0("bray_dist ~ ", group_col))
# }
# 
# dbRDA_model <- capscale(dbRDA_formula, data = meta)
# saveRDS(dbRDA_model, file.path(out_dir, "dbrda_model.rds"))
# 
# anova_all <- anova(dbRDA_model, permutations = cfg$differential$dbrda$permutations)
# anova_margin <- anova(dbRDA_model, by = "margin", permutations = cfg$differential$dbrda$permutations)
# 
# write_csv_safe(as.data.frame(anova_all) %>% tibble::rownames_to_column("Term"), file.path(out_dir, "dbrda_anova_all.csv"))
# write_csv_safe(as.data.frame(anova_margin) %>% tibble::rownames_to_column("Term"), file.path(out_dir, "dbrda_anova_margin.csv"))
# 
# # CHANGE C: formato limpio y parseable de R²
# r2 <- RsquareAdj(dbRDA_model)
# writeLines(c(
#   paste0("R2:           ", r2$r.squared),
#   paste0("Adjusted R2:  ", r2$adj.r.squared)
# ), con = file.path(out_dir, "dbrda_r2.txt"))
# 
# # envfit: prevalence filter
# otu_mat <- as(otu_table(pseq), "matrix")
# if (!taxa_are_rows(pseq)) otu_mat <- t(otu_mat)
# otu_mat <- t(otu_mat)[common_samples, , drop = FALSE]
# 
# prev <- apply(otu_mat, 2, function(x) sum(x > 0) / length(x))
# otu_f <- otu_mat[, prev > cfg$differential$dbrda$prevalence_min, drop = FALSE]
# 
# # CHANGE B: p.max = 1 — crítico para que pvals se alineen con scores(fit).
# # Sin esto, envfit descarta vectores no significativos y pvals quedan
# # reciclados/misalineados → falsos positivos en envfit_significant_otus.csv.
# fit <- envfit(dbRDA_model, otu_f,
#               perm  = cfg$differential$dbrda$permutations,
#               p.max = 1)
# 
# fit_scores <- as.data.frame(scores(fit, display = "vectors"))
# fit_scores$pval <- fit$vectors$pvals
# fit_scores$OTU_ID <- rownames(fit_scores)
# 
# # Sanity check: scores y pvals deben tener el mismo largo
# stopifnot(nrow(fit_scores) == length(fit$vectors$pvals))
# 
# fit_sig <- fit_scores[!is.na(fit_scores$pval) & fit_scores$pval < 0.05, , drop = FALSE]
# 
# tax_df <- as.data.frame(tax_table(pseq)) %>% tibble::rownames_to_column("OTU_ID")
# fit_sig_annot <- merge(fit_sig, tax_df, by = "OTU_ID")
# 
# write_csv_safe(fit_sig, file.path(out_dir, "envfit_significant_otus.csv"))
# write_csv_safe(fit_sig_annot, file.path(out_dir, "envfit_significant_otus_with_taxonomy.csv"))
# 
# # 기본 플롯 저장(PDF 권장)
# pdf(file.path(out_dir, "dbrda_sites.pdf"), width = 6, height = 5)
# plot(dbRDA_model, display = "sites", main = "db-RDA")
# dev.off()
# 
# cat("✅ Script 07 finalizado exactamente según el diseño original.\n")
# 


