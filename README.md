#READ.ME 12NOV26

# 🧬 16s_insomnia: Microbiota Analysis in Insomnia

Bienvenido al repositorio oficial del pipeline bioinformático y de metaanálisis para el estudio del rol de la microbiota intestinal en el insomnio. Este proyecto replica y expande el diseño metodológico de Park et al. (2026), integrando un flujo robusto que abarca desde el preprocesamiento de lecturas crudas de secuenciación de amplicomes 16S rRNA hasta análisis avanzados de diversidad, aprendizaje automático (Boruta), modelado multivariado (MaAsLin2, db-RDA, LEfSe) y predicción funcional mediante PICRUSt2.

---

## 🏗️ Estructura General del Proyecto

```text
16s_insomnia/
├── 16s_insomnia.Rproj
├── config/
│   └── config.yaml                               <-- Archivo de configuración central[cite: 5]
│
├── 0_preprocessing/
│   ├── preprocessing_README.md
│   ├── preprocessing_README.html
│   ├── config/
│   │   └── datasets_config.tsv                   <-- Dynamic parameters (primers, truncation, layout)
│   ├── notes/
│   ├── script_drafts_ideas/                      <-- Drafts and experimental tests
│   └── scripts/
│       ├── 00_install_sra_tools.sh               <-- SRA tools installation
│       ├── 01_download_silva_classifier.sh       <-- SILVA classifier download (Full-Length)
│       ├── 01_download_study_sras_w_metadata.sh  <-- Raw FASTQ and SRA metadata download
│       ├── 01_install_qiime2.sh                  <-- QIIME 2 environment setup
│       ├── 02_metadata_qiime2_format.sh          <-- Base metadata formatting for QIIME 2
│       ├── 03_merge_study_w_downloaded_metadata.py <-- Enrich and harmonize clinical metadata (ISI, BMI, etc.)
│       ├── 04_qiime_import_and_trim.sh           <-- FASTQ import and primer removal (Cutadapt)
│       ├── 04_viz_quality_plot.sh                <-- .qzv generation for quality inspection
│       ├── 05_dada2_denoise_datasets.sh          <-- Denoising and error correction with DADA2
│       ├── 05b_visualize_denoise_stats.sh        <-- Statistical summary of DADA2 process
│       ├── 06_merge_cluster_export.sh            <-- Merge, VSEARCH 99% clustering, SILVA taxonomy, and R export
│       ├── 06b_visualize_unified_table.sh        <-- Generates summary visualization of unified table (.qzv)
│       └── 06D_build_master_metadata.py          <-- Builds final master metadata table (master_metadata_unified.tsv)
│
├── utils/                                        <-- Funciones auxiliares de R (raíz)
│   ├── utils_io.R
│   ├── utils_phyloseq.R
│   ├── utils_plot.R
│   └── utils_stats.R
│
├── 1_meta_analysis/                              <-- Scripts del pipeline analítico (01 a 15)
│
├── data/
│   ├── raw/                                      <-- Raw FASTQs per study and supplementary tables
│   └── pre_processed/
│       ├── metadata/                             <-- Formatted (*_q2.tsv) and enriched metadata
│       ├── unified/                              <-- Inputs for meta-analysis: unified_otu_table.tsv, 
│       │                                             taxonomy.tsv, master_metadata_unified.tsv, 
│       │                                             unified_rep_seqs.fasta
│       └── qiime2_workspace/
│
├── resources/                                    <-- SILVA classifiers and external resources
│   └── classifiers/
│
├── picrust2_out/                                 <-- Carpeta de salida temporal/interna de PICRUSt2
│
└── outputs/                                      <-- ¡Resultados oficiales del metaanálisis!
    ├── tables/
    ├── models/
    ├── logs/
    ├── figures/
    ├── dbrda/
    ├── lefse/
    ├── boruta/
    └── function/

```

---

## ⚙️ Entorno de Trabajo y Prerrequisitos

Todo el pipeline se ejecuta bajo un entorno virtual centralizado de **Conda** optimizado para bioinformática (`qiime2-amplicon-2024.5`).

### Instalación Rápida con el Script Maestro (`setup_and_run.sh`)

Para desplegar y configurar automáticamente el proyecto tanto en tu PC con **Linux Mint** como en **macOS**, utiliza el asistente automatizado ubicado en la raíz:

1. Abre tu terminal en la carpeta principal del proyecto (`16s_insomnia/`).
2. Otorga permisos de ejecución (solo la primera vez):

```bash
chmod +x setup_and_run.sh
```

3. Ejecuta el script maestro:

```bash
./setup_and_run.sh
```

Este script detectará tu sistema, validará el entorno Conda `qiime2-amplicon-2024.5` y sincronizará de forma transparente todas las dependencias de R utilizando **`renv`** (`renv.lock`).

---

## 📥 0. Módulo de Preprocesamiento (FASTQ → Tablas Unificadas)

> **Nota importante sobre el flujo:** Este pipeline está diseñado para integrar inicialmente **3 estudios principales** (con la capacidad prevista de incorporar otros **5 estudios adicionales** de forma modular). Debido a esta integración progresiva y a la heterogeneidad de los metadatos públicos, el preprocesamiento no cuenta con un único script de ejecución global (`run_all`), requiriendo una supervisión y ajuste manual de los metadatos por cohorte.

### 📋 Resumen del Flujo de Preprocesamiento

1. **Adquisición:** Descarga de archivos FASTQ crudos y metadatos desde SRA utilizando SRA Toolkit.

2. **Armonización:** Normalización de variables clínicas clave (escalas ISI, IMC, regiones geográficas y estandarización de los grupos de estudio en `Insomnia` vs `Control`).

3. **Procesamiento QIIME 2 / DADA2:**
   - Remoción de primers específicos de región mediante *Cutadapt*.
   - Corrección de errores y *denoising* con *DADA2* (ajustando dinámicamente parámetros de truncamiento *Paired-End* o *Single-End*).
   - Filtrado de muestras de baja profundidad (< 10,000 lecturas retenidas).
   - Clasificación taxonómica mediante un clasificador Naive Bayes entrenado en la base de referencia **SILVA 138.1** (identidad al 99%).

4. **Integración Global:** Fusión de tablas de características mediante *VSEARCH* (*de novo* clustering al 99% de identidad) para generar un vocabulario unificado de ASVs/OTUs compartido entre todas las cohortes, exportando matrices limpias listas para el análisis en R.

Para una revisión detallada de cada script, parámetros de configuración (`datasets_config.tsv`) y diagramas de ejecución de esta fase, consulta la documentación específica en [0_preprocessing/README.md](0_preprocessing/README.md) o revisa la guía interna de la cohorte.

---

## 📊 1. Módulo de Metaanálisis y Análisis Funcional en R (`1_meta_analysis/`)

Una vez consolidada la matriz unificada (`unified_otu_table.tsv`) y el maestro de metadatos (`master_metadata_unified.tsv`), el flujo estadístico avanzado se ejecuta de manera secuencial mediante el script maestro:

```bash
conda activate qiime2-amplicon-2024.5
Rscript 1_meta_analysis/run_all.R
```

El pipeline de metaanálisis ejecuta de forma automatizada los siguientes componentes (respaldados por el diseño optimizado de Park et al.):

- **01. Construcción de Objeto Phyloseq (`01_build_phyloseq.R`):** Carga las tablas unificadas, filtra secuencias no deseadas (mitocondrias, cloroplastos o asignaciones nulas) y fija los niveles de los factores de estudio.
- **02. Diversidad Alfa (`02_alpha_diversity.R`):** Calcula índices de riqueza y diversidad (Shannon, Simpson) utilizando gráficos de violín/caja compuestos mediante *patchwork* y pruebas estadísticas t-test.
- **03. Diversidad Beta y PERMANOVA (`03_beta_diversity_permanova.R`):** Evalúa la disimilitud de Bray-Curtis ajustando por covariables clínicas multivariadas (edad, sexo, cohorte/estudio) mediante *vegan::adonis2*.
- **04. Razón Firmicutes / Bacteroidota (`04_fb_ratio.R`):** Cuantifica de forma robusta la proporción de los principales filos microbianos asociados al sueño.
- **05 & 06. Selección de Características por Machine Learning - Boruta (`05_boruta.R` & `06_boruta_pca_heatmap.R`):** Identifica ASVs discriminantes clave e integra análisis de componentes principales (PCA) junto a mapas de calor normalizados por Z-score.
- **07. Ordenación Restringida db-RDA y Envfit (`07_dbrda_envfit_fix.R`):** Ejecuta selección de modelos escalonados (*ordiR2step*) y correlaciona la matriz comunitaria con variables ambientales/clínicas corrigiendo desalineaciones de p-valores.
- **08. Biomezcladores diferenciales LEfSe (`08_lefse.R`):** Detecta biomarcadores taxonómicos robustos utilizando la API moderna de *microbiomeMarker* y genera gráficos de barras con puntuaciones LDA firmadas.
- **09. Abundancia Diferencial MaAsLin2 (`09_maaslin2_taxa.R`):** Aplica modelos lineales multivariados evitando colapsos gráficos en versiones modernas de ggplot2 mediante un bypass visual explícito.
- **10 al 15. Predicción Funcional (PICRUSt2 y Redes Metabólicas):**
  - Convierte y procesa perfiles metagenómicos predichos hacia ortólogos KEGG (KOs).
  - Calcula correlaciones vectorizadas de Spearman entre taxones clave y funciones metabólicas.
  - Genera diagramas de acorde (*Chord Diagrams* optimizados con límites de memoria RAM) y diagramas de puntos (*Dotplots*) para la visualización de rutas enriquecidas.

---

## 📌 Puntos Clave para Reportes y Tesis

1. **Robustez en Preprocesamiento:** El uso de una base unificada por VSEARCH (99%) previene sesgos derivados de llamadas de ASVs independientes por estudio.

2. **Control de Versiones y Compatibilidad:** La desactivación de motores gráficos internos obsoletos en MaAsLin2 y el manejo defensivo de memoria en diagramas circulares garantizan una ejecución libre de fallos OOM (Out Of Memory) en equipos locales.

3. **Reproducibilidad Total:** Gracias al versionamiento de entornos virtuales y la gestión automatizada de librerías mediante *renv*, cualquier análisis puede replicarse de forma idéntica en distintas arquitecturas (Linux Mint / macOS).

---

## Installation & Environment Setup

This project integrates QIIME 2 and the R package `phyloseq`. To prevent environmental dependency conflicts inside the Conda environment, R dependencies must be built natively using system headers.

### 1. Prerequisites (Linux System Libraries)
Before setting up the environment, you must install the required development headers for compiling data transfer and XML manipulation packages in Linux:
```bash
sudo apt update
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev
```

### 2. Activate Environment
Activate your managed QIIME 2 deployment containing your pre-compiled dependencies:
```bash
conda activate qiime2-amplicon-2024.5
```

### 3. Install Phyloseq (Inside the Environment)
To keep the environment stable, launch the internal R environment directly from your terminal session and deploy `phyloseq` using `BiocManager`:

1. Open the environment-bound R console:
   ```bash
   R
   ```
2. Run the deployment execution block within the active session:
   ```R
   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("phyloseq")
   q(save="no")
   ```

### 4. Verification
Verify that the package initializes successfully without throwing compilation errors:
```bash
Rscript -e "library(phyloseq); cat('✅ Phyloseq initialization verified!\n')"
```
---

---
# Pipeline: Inputs → Outputs per Script

## 1_meta_analysis/
# Pipeline: Inputs → Outputs per Script (`1_meta_analysis/`)

---

### `00_setup.R`

* **Inputs:**
* `config/config.yaml` (from user configuration setup)


* **Outputs:** Loaded global environment configuration (`cfg`), random seed assignment, and shared R library environments.

---

### `01_build_phyloseq.R`

* **Inputs:**
* `data/pre_processed/unified/unified_otu_table.tsv` (from `06_merge_cluster_export.sh`)
* `data/pre_processed/unified/taxonomy.tsv` (from `06_merge_cluster_export.sh`)
* `data/pre_processed/unified/master_metadata_unified.tsv` (from `06D_build_master_metadata.py`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/models/physeq_filtered.rds` ⚠️ **Primary filtered phyloseq object**
* `outputs/tables/otu_filtered.csv`
* `outputs/tables/tax_filtered.csv`
* `outputs/tables/metadata_filtered.csv`



---

### `02_alpha_diversity.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/tables/alpha_diversity.csv`
* `outputs/figures/alpha_shannon_simpson.png`
* `outputs/logs/alpha_diversity_stats.txt`



---

### `03_beta_diversity_permanova.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/tables/beta_diversity.csv`
* `outputs/figures/beta_pcoa_nmds.png`
* `outputs/logs/permanova.txt`



---

### `04_fb_ratio.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/tables/fb_ratio.csv`
* `outputs/figures/fb_ratio.png`
* `outputs/logs/fb_ratio_stats.txt`



---

### `05_boruta.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/boruta/boruta_model.rds`
* `outputs/boruta/boruta_selected_otu.csv`
* `outputs/boruta/boruta_selected_otus_with_taxonomy.csv`
* `outputs/boruta/boruta_importance_clean.png`



---

### `06_boruta_pca_heatmap.R` (+ QC script)

* **Inputs:**
* `outputs/boruta/boruta_selected_otu.csv` (from `05_boruta.R`)
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/boruta/boruta_pca.png`
* `outputs/boruta/boruta_heatmap_samples.png`



---

### `07_dbrda_envfit_fix.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/dbrda/ordiR2step_model.rds`
* `outputs/dbrda/dbrda_model.rds`
* `outputs/dbrda/dbrda_anova_all.csv`
* `outputs/dbrda/dbrda_anova_margin.csv`
* `outputs/dbrda/dbrda_r2.txt`
* `outputs/dbrda/envfit_significant_otus_with_taxonomy.csv`
* `outputs/dbrda/dbrda_sites.pdf`



---

### `08_lefse.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/lefse/lefse_result.csv`
* `outputs/lefse/lefse_sig_with_otu.csv`
* `outputs/lefse/lefse_barplot.png`



---

### `09_maaslin2_taxa.R`

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/maaslin2/` (MaAsLin2 output directory)
* `outputs/maaslin2/maaslin2_sig_with_taxonomy.csv`



---

### `10A_run_picrust2.sh`

* **Inputs:**
* `data/pre_processed/unified/unified_rep_seqs.fasta` (from `06_merge_cluster_export.sh`)
* `data/pre_processed/unified/unified_otu_table.tsv` (or `.biom`) (from `06_merge_cluster_export.sh`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/tables/merged_KO_metagenome.tsv` ⚠️ **Predicted unstratified functional profiles**



---

### `10B_function_taxa_ko_cor.R` (+ QC script)

* **Inputs:**
* `outputs/tables/merged_KO_metagenome.tsv` (from `10A_run_picrust2.sh`)
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/function/tables/taxon_ko_cor_all.csv`
* `outputs/function/tables/taxon_ko_edges_filtered.csv` ⚠️ **Significant Spearman taxon-KO correlations**



---

### `11_function_chord_edges.R`

* **Inputs:**
* `outputs/function/tables/taxon_ko_edges_filtered.csv` (from `10B_function_taxa_ko_cor.R`)
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/function/figures/chord_taxon_ko.png`
* `outputs/function/figures/chord_taxon_ko.pdf`



---

### `12_function_dotplot.R`

* **Inputs:**
* `outputs/function/tables/taxon_ko_edges_filtered.csv` (from `10B_function_taxa_ko_cor.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/function/figures/figure5b_taxon_ko_dotplot.png`



---

### `13a_extract_ko_contrib.R` (+ QC script)

* **Inputs:**
* `outputs/models/physeq_filtered.rds` (from `01_build_phyloseq.R`)
* `outputs/function/tables/taxon_ko_edges_filtered.csv` (from `10B_function_taxa_ko_cor.R`)
* `outputs/tables/merged_KO_metagenome.tsv` (from `10A_run_picrust2.sh`)


* **Outputs:**
* `outputs/tables/metadata_for_taxa.csv`
* `outputs/tables/picrust_taxonomy_filled.csv`
* `outputs/tables/KO_FC_group.csv`
* `outputs/tables/ko_contrib_filtered.csv` ⚠️ **Filtered per-sample functional contributions**



---

### `13b_function_contribution_diff.R` (+ QC script)

* **Inputs:**
* `outputs/tables/ko_contrib_filtered.csv` (from `13a_extract_ko_contrib.R`)
* `outputs/tables/metadata_for_taxa.csv` (from `13a_extract_ko_contrib.R`)
* `outputs/tables/KO_FC_group.csv` (from `13a_extract_ko_contrib.R`)


* **Outputs:**
* `outputs/function_contrib/differential_functional_contributions.csv`



---

### `14_export_ko_list_for_cluego.R`

* **Inputs:**
* `outputs/function/tables/taxon_ko_edges_filtered.csv` (from `10B_function_taxa_ko_cor.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/function_contrib/network/cluego_gene_list_KO.txt`



---

### `15_function_chord_from_ko_fc.R`

* **Inputs:**
* `outputs/tables/KO_FC_group.csv` (from `13a_extract_ko_contrib.R`)
* `config/config.yaml` (from user configuration setup)


* **Outputs:**
* `outputs/function/figures/chord_ko_fc.png`



---

### `run_all.R`

* **Inputs:** Sequential execution queue of all active preprocessing and metaanalysis scripts (`01` through `15`).
* **Outputs:** Full runtime logs and complete unified `outputs/` directory structure.
---

---
## 0_preprocessing/

---

### `00_install_sra_tools.sh`
- **Inputs:** None (checks system for `prefetch` / `fasterq-dump`)
- **Outputs:** SRA Toolkit installed at `/usr/local/sra-toolkit/bin` (system-level; no repo files)

---

### `01_install_qiime2.sh`
- **Inputs:** None (downloads `qiime2-amplicon-2024.5_py39-osx-conda.yml` from data.qiime2.org)
- **Outputs:** Conda environment `qiime2-amplicon-2024.5`

---

### `01_download_silva_classifier.sh`
- **Inputs:** Official URL `silva-138-99-nb-classifier.qza` (~500 MB)
- **Outputs:**
  - `resources/classifiers/silva-138-99-nb-classifier.qza` ⚠️ **critical resource**

---

### `01b_download_study_sras.sh`
- **Inputs:** `0_preprocessing/config/sra_projects.tsv` (DATASET, PROJECT columns)
- **Outputs:**
  - `data/raw/{run_id}/*.fastq.gz` (compressed FASTQs)
  - SRA metadata downloaded via `pysradb`

---

### `02_metadata_qiime2_format.sh`
- **Inputs:**
  - `0_preprocessing/config/datasets_config.tsv`
  - `data/raw/{STUDY}/{STUDY}_metadata_detailed.csv` (from step 01b)
  - `0_preprocessing/utils/02_metadata_qiime2_format_engine.py`
- **Outputs:**
  - `data/preprocessed/metadata/{STUDY}_metadata_q2.tsv` ⚠️ **QIIME 2-formatted metadata (with `#q2:types` header)**

---

### `03_merge_study_w_downloaded_metadata.py`
- **Inputs:**
  - `data/preprocessed/metadata/{STUDY}_metadata_q2.tsv`
  - `data/raw/{STUDY}/table_s1.csv` (supplementary clinical data, study-dependent)
  - `0_preprocessing/utils/metadata_q2_enriched_engine.py`
- **Outputs:**
  - `data/preprocessed/metadata/{STUDY}_metadata_q2_format_enriched.tsv` ⚠️ **enriched metadata ready for merging**

---

### `04_qiime_import_and_trim.sh`
- **Inputs:**
  - `0_preprocessing/config/datasets_config.tsv` (primers, truncation params)
  - `data/raw/{STUDY}/*.fastq.gz`
- **Outputs:**
  - `data/pre_processed/qiime2_workspace/manifests/{STUDY}_manifest.tsv`
  - `data/pre_processed/qiime2_workspace/{STUDY}/demux.qza`
  - `data/pre_processed/qiime2_workspace/{STUDY}/demux_trimmed.qza` ⚠️ **primer-trimmed reads**

---

### `04_viz_quality_plot.sh`
- **Inputs:**
  - `data/pre_processed/qiime2_workspace/{STUDY}/demux_trimmed.qza`
- **Outputs:**
  - `data/pre_processed/qiime2_workspace/{STUDY}/demux_trimmed.qzv` (quality visualization)

---

### `05_dada2_denoise_datasets.sh`
- **Inputs:**
  - `0_preprocessing/config/datasets_config.tsv` (truncation lengths, paired/single flag)
  - `data/pre_processed/qiime2_workspace/{STUDY}/demux_trimmed.qza`
- **Outputs (per study):**
  - `{STUDY}/table.qza` ⚠️ **feature table**
  - `{STUDY}/rep_seqs.qza` ⚠️ **representative sequences**
  - `{STUDY}/denoising_stats.qza`

---

### `05b_visualize_denoise_stats.sh`
- **Inputs:**
  - `0_preprocessing/config/datasets_config.tsv`
  - `data/pre_processed/metadata/{STUDY}_metadata_q2.tsv`
  - `{STUDY}/denoising_stats.qza`, `{STUDY}/table.qza`, `{STUDY}/rep_seqs.qza`
- **Outputs:**
  - `{STUDY}/denoising_stats.qzv`
  - `{STUDY}/table.qzv`
  - `{STUDY}/rep_seqs.qzv`

---

### `06_merge_cluster_export.sh`
- **Inputs:**
  - `0_preprocessing/config/datasets_config.tsv`
  - `resources/classifiers/silva-138-99-nb-classifier.qza`
  - `qiime2_workspace/{STUDY}/table.qza` + `rep_seqs.qza` (all active studies)
- **Outputs:**
  - **Internal (QIIME 2 workspace):**
    - `qiime2_workspace/unified/merged_table_raw.qza`
    - `qiime2_workspace/unified/merged_rep_seqs.qza`
    - `qiime2_workspace/unified/final_unified_table.qza` ⚠️ **VSEARCH 99% clustered**
    - `qiime2_workspace/unified/unified_rep_seqs.qza`
    - `qiime2_workspace/unified/taxonomy.qza`
  - **Human-readable export:**
    - `data/pre_processed/unified/unified_otu_table.tsv` ⚠️
    - `data/pre_processed/unified/unified_otu_table.csv` ⚠️
    - `data/pre_processed/unified/unified_rep_seqs.fasta` ⚠️
    - `data/pre_processed/unified/taxonomy.tsv` ⚠️

---

### `06b_visualize_unified_table.sh`
- **Inputs:**
  - `qiime2_workspace/unified/final_unified_table.qza`
  - `qiime2_workspace/unified/unified_rep_seqs.qza`
  - `qiime2_workspace/unified/taxonomy.qza`
- **Outputs:**
  - `qiime2_workspace/unified/final_unified_table.qzv`
  - `qiime2_workspace/unified/unified_rep_seqs.qzv`
  - `qiime2_workspace/unified/taxonomy.qzv`

---

### `06D_build_master_metadata.py`
- **Inputs:**
  - `data/preprocessed/metadata/*_metadata_q2_format_enriched.tsv` (glob)
- **Outputs:**
  - `data/preprocessed/unified/master_metadata_unified.tsv` ⚠️ **master metadata for R/phyloseq**

---

## Key "hand-off" artifacts to `1_meta_analysis/`

| Artifact | Produced by | Consumed by |
|---|---|---|
| `unified_otu_table.tsv` | `06_merge_cluster_export.sh` | `utils_io.R` → `01_build_phyloseq.R` |
| `taxonomy.tsv` | `06_merge_cluster_export.sh` | `utils_io.R` → `01_build_phyloseq.R` |
| `unified_rep_seqs.fasta` | `06_merge_cluster_export.sh` | `10A_run_picrust2.sh` |
| `master_metadata_unified.tsv` | `06D_build_master_metadata.py` | `utils_io.R` → `01_build_phyloseq.R` |
| `merged_KO_metagenome.tsv` | `10A_run_picrust2.sh` | `13a_extract_ko_contrib.R` |


