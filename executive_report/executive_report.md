# Executive Report: 16S rRNA Microbiome Meta-Analysis in Insomnia Disorder

**Repository:** `https://github.com/Dskies5i/16s_insomnia`  
**Date:** November 2026  
**Cohort Scope (Initial Meta-Analysis Phase):** $N = 139$ samples (**Control:** $n = 25$; **Insomnia:** $n = 114$) across three harmonized studies (`chen_2022`, `haimov_2022`, `liub_2019`), consolidated into **2,394 filtered ASVs/taxa** (VSEARCH 99% clustering, SILVA 138 classification, after removing `g_other`, `Unassigned`, and zero-sum features).

---

## 1. Executive Summary & Pipeline Status

* **End-to-End Pipeline Validation:** Scripts `00` through `15` (including DADA2 denoising, VSEARCH 99% global clustering, SILVA classification, phyloseq filtering, diversity metrics, Boruta feature selection, constrained `db-RDA`/`envfit`, LEfSe, MaAsLin2, and PICRUSt2 functional inference + Spearman taxon–KO networks) have been validated end-to-end on macOS, with the previously blocking `10A_run_picrust2.sh` step resolved.
* **Preprocessing Added Beyond Original Pipeline:** The **original pipeline had no preprocessing module**; `0_preprocessing/` was added to standardize SRA metadata (`#q2:types`) and clinical covariates, paired with `setup_and_run.sh` and `renv` for cross-platform execution across macOS and Linux environments.
* **Core Biological Signal:** Insomnia samples show reduced alpha diversity, a lower *Firmicutes*:*Bacteroidota* ratio, a pronounced shift toward depletion of commensal taxa in LEfSe (44/56 biomarkers enriched in Control vs. 12/56 in Insomnia), 21 Boruta-selected discriminant ASVs, and broad predicted KEGG Ortholog (KO) remodeling (3,530 upregulated KOs in Insomnia vs. 5,412 downregulated/Control-preserved KOs).

---

## 2. Pipeline Architecture & Engineering Walkthrough

### Dataset & Cohort Harmonization
* **Multi-Cohort 16S rRNA Workflow (Control vs. Insomnia):** Harmonizes sequencing data and clinical phenotypes across multiple studies (`chen_2022`, `haimov_2022`, `liub_2019`) into a unified Control vs. Insomnia comparison framework.
* **QIIME 2 / DADA2 Denoising & VSEARCH 99% Clustering:** Processes raw amplicon reads through DADA2 denoising and VSEARCH 99% clustering with SILVA 138 taxonomic assignment and phyloseq filtering (`g_other`, `Unassigned`, and zero-sum features removed).
* **Clinical Covariate Adjustment (`study_id`, `age`, `sex`):** Accounts for cohort/batch structure (`study_id`) alongside host covariates (`age`, `sex`) and clinical metadata so downstream statistical and ordination models separate disease effects from study or demographic confounding.

### Analytical & ML Stack
* **Alpha/Beta Diversity + PERMANOVA:** Evaluates within-sample richness/evenness (including Shannon and Simpson) and between-sample community dissimilarity (Bray–Curtis) paired with PERMANOVA and ordination visualizations.
* ***Firmicutes*:*Bacteroidota* Ratio:** Computes phylum-level *Firmicutes*:*Bacteroidota* ratios across samples and tests group differences between Control and Insomnia cohorts.
* **Boruta Feature Selection + PCA/Z-Score Heatmaps:** Uses Boruta random-forest feature selection to identify discriminant ASVs and visualizes sample/feature structure via PCA and Z-score heatmaps.
* **Constrained `db-RDA` + `envfit`:** Models constrained distance-based redundancy analysis (`db-RDA`) and fits environmental/clinical vectors (`envfit`) to evaluate how covariates and clinical variables relate to community composition.
* **LEfSe & MaAsLin2 (with Modern `ggplot2` Plotting Safeguards):** Combines LEfSe biomarker discovery with covariate-adjusted MaAsLin2 differential abundance testing, updated with modern `ggplot2` plotting safeguards to prevent breakage from package updates.
* **PICRUSt2 Functional/KO Correlation & Contribution Analysis:** Infers metagenomic functional potential (KEGG Orthologs) via PICRUSt2, quantifies group fold changes, computes taxon–KO correlations, and extracts functional contribution tables (`13a`/`13b`) and CluEgo exports (`14`).

### Engineering & Reproducibility Fixes
* **Vectorized Spearman Correlations:** Replaces slow loop-based correlation calculations with vectorized computation, cutting runtime from hours to seconds.
* **Memory-Safe Chord Diagram Capping:** Caps edge counts prior to rendering chord diagrams to prevent out-of-memory (OOM) crashes on large taxon–KO correlation tables.
* **Cross-Platform `setup_and_run.sh` + `renv` Setup:** Adds a unified `setup_and_run.sh` workflow and `renv` lockfile/environment management to support reproducible runs across macOS and Linux.

---

## 3. Biological Interpretation

### A. Loss of Ecological Diversity & Commensal Depletion
* **Reduced Alpha Diversity (Shannon & Simpson):** Insomnia is associated with a statistically significant decline in both Shannon diversity ($p = 0.016$; mean $3.253$ in Insomnia vs. $3.530$ in Control; Welch's $t = 2.51$) and Simpson diversity ($p = 0.023$; mean $0.907$ in Insomnia vs. $0.931$ in Control; Welch's $t = 2.32$). In gut ecology, lower richness and evenness reflect a less resilient microbial community that is more vulnerable to pathobiont expansion and metabolic perturbation.
* **Asymmetric Biomarker Shift Toward Commensal Loss:** LEfSe identified **56 significant genus-level biomarkers** ($\text{LDA} \ge 2.037$, mean $\text{LDA} = 2.994$, $p < 0.05$), of which **44 taxa are enriched in Controls** and only **12 in Insomnia**. Rather than being defined solely by a single "insomnia pathogen," the insomnia microbiome exhibits a broad contraction of baseline commensal diversity—a pattern frequently linked in gut–brain axis research to altered neuroimmune tone and mucosal homeostasis.

### B. *Firmicutes*:*Bacteroidota* Ratio & Key Taxonomic Signatures
* **Shift in *Firmicutes*:*Bacteroidota* ($N = 138$ valid samples):** Welch's $t$-test shows a significant difference between groups ($t = 2.77$, $p = 0.01$), with a lower mean ratio in Insomnia ($1.142$) compared to Control ($1.952$) and a cohort median of $1.096$. Because broad phylum ratios aggregate ecologically diverse members, this shift is best interpreted alongside fine-grained ASV and genus-level signals.
* **SCFA-Producing & Gut–Brain Axis Taxa (Boruta & LEfSe):** Machine learning feature selection via Boruta confirmed **21 discriminant ASVs**, highlighting members of *Lachnospiraceae* and *Erysipelotrichaceae* alongside genera such as *Coprococcus*, *Roseburia*, *Butyricicoccus*, *[Eubacterium] hallii group*, *Dorea*, *Agathobacter*, and *Bacteroides*. Many *Lachnospiraceae* and butyrate-producing genera support intestinal barrier integrity, histone deacetylase (HDAC) inhibition, and vagal/neuroinflammatory regulation; reductions or compositional shifts in short-chain fatty acid (SCFA) producers are plausible contributors to sleep fragmentation and low-grade systemic inflammation.
* **Potential Opportunistic / Pathobiont & Oral/Upper-GI Signals:** The presence of genera such as *Streptococcus*, *Actinomyces*, *Gemella*, *Erysipelatoclostridium*, or *Acinetobacter* among candidate/discriminant features warrants careful evaluation against cohort and sequencing-region effects, as shifts in non-commensal or oral-associated taxa in stool profiles can reflect altered gut transit, mucosal stress, or study-specific batch structure.

### C. Community Ordination, Batch Effects & Covariate Confounding
* **Bray–Curtis Ordination & Constrained Modeling (`db-RDA`):** PCoA and NMDS show broad visual overlap with directional dispersion in Insomnia samples, while the conditioned `db-RDA` model adjusting for clinical/technical covariates yields a low explained variance ($R^2 = 0.0062$, adjusted $R^2 = 0.0003$).
* **Importance of Controlling `study_id` and Technical Heterogeneity:** Both MaAsLin2 and `db-RDA` underscore strong cohort/batch structure (`study_id`, including `haimov_2022` and differences across V3–V4 vs. V4 amplicon regions and single- vs. paired-end truncation lengths). In multi-cohort 16S meta-analyses, unadjusted group comparisons can conflate study origin, primer bias, and clinical covariates (such as age, sex, BMI, or insomnia severity scales like ISI) with true disease biology—making covariate-adjusted models (`study_id`, `age`, `sex`) essential for separating genuine host–microbe signals from technical artifacts.

### D. Predicted Functional Remodeling (PICRUSt2 & Taxon–KO Networks)
* **Asymmetric KO Regulation:** Across **8,942 predicted KEGG Orthologs (KOs)**, **3,530 KOs are upregulated in Insomnia** ($\text{FC} > 0$) whereas **5,412 KOs are downregulated / higher in Control** ($\text{FC} < 0$). This mirrors the taxonomic contraction in Insomnia: a broad loss of functional capacity preserved in healthy controls alongside selective upregulation of metabolic functions in insomnia.
* **Dense Taxon–Function Coupling:** Evaluating **4,533,594 taxon–KO pairs** (507 unique taxa $\times$ 8,942 KOs) yielded **76,537 significant Spearman correlations** ($|\rho| \ge 0.5$, $1.69\%$ of pairs; $\rho$ range $-0.692$ to $1.000$, mean $|\rho| = 0.665$). Coupled with differential contribution tables (`13a`/`13b`) and CluEgo KO exports (`14`), these associations provide candidates for pathway-level enrichment (e.g., amino acid/neurotransmitter precursors, SCFA and bile acid metabolism, and stress/oxidative pathways) while keeping in mind that PICRUSt2 infers gene content from 16S phylogeny rather than direct metagenomic sequencing.

---

## 4. Summary of Key Quantitative Findings

| Metric / Module | Key Finding | Notes & Context |
| :--- | :--- | :--- |
| **Sample & Feature Scope** | $N = 139$ (Control $25$, Insomnia $114$); **2,394 ASVs** | F:B ratio evaluated on $N = 138$ (1 sample lacked reads in target phyla) |
| **Alpha Diversity** | Shannon $p = 0.016$ ($3.530$ vs. $3.253$); Simpson $p = 0.023$ ($0.931$ vs. $0.907$) | Lower diversity in Insomnia via Welch's two-sample $t$-test |
| **F:B Ratio** | $p = 0.01$ ($t = 2.77$; mean $1.952$ Control vs. $1.142$ Insomnia) | Cohort median $1.096$ (range $0.133\text{–}6.401$) |
| **Boruta ML** | **21 confirmed discriminant ASVs** | Enriched for *Lachnospiraceae* and *Erysipelotrichaceae* members |
| **LEfSe** | **56 significant biomarkers** ($\text{LDA} \ge 2.037$, mean $\text{LDA} = 2.994$, $p < 0.05$) | **44** enriched in Control, **12** enriched in Insomnia |
| **PICRUSt2 & Spearman** | **8,942 KOs** (**3,530** $\text{FC} > 0$, **5,412** $\text{FC} < 0$); **76,537** significant edges ($|\rho| \ge 0.5$) | Out of $4,533,594$ tested taxon–KO pairs (507 taxa $\times$ 8,942 KOs) |

---

## 5. Pipeline Scripts, Preprocessing & Key Improvements

*Note: The original pipeline had **no preprocessing module** (`0_preprocessing/`), and required environment/runtime hardening for multi-cohort reproducibility and large functional correlation matrices.*

| Pipeline Component / Stage | Role in Workflow | Key Additions & Improvements Made |
| :--- | :--- | :--- |
| **`0_preprocessing/`** | Metadata & covariate standardization | **New module** (original pipeline had no preprocessing); standardizes SRA metadata (`#q2:types`) and clinical covariates (`study_id`, `age`, `sex`) |
| **`setup_and_run.sh` & `renv`** | Environment & execution orchestration | Added cross-platform setup and locked R package versions (`renv`) for reproducible runs across macOS and Linux |
| **Scripts `00`–`09` (QIIME 2 / DADA2, VSEARCH, SILVA, Diversity, Boruta, `db-RDA`, LEfSe, MaAsLin2)** | Amplicon processing, diversity, ML feature selection & differential abundance | Validated end-to-end on macOS; incorporated covariate adjustments (`study_id`, `age`, `sex`) and modern `ggplot2` plotting safeguards |
| **`10A_run_picrust2.sh` & Scripts `10`–`15` (PICRUSt2, KO fold changes, Spearman networks, `13a`/`13b`, CluEgo `14`)** | Functional inference, taxon–KO correlations & contribution summaries | Resolved blocking `10A_run_picrust2.sh` step; vectorized Spearman correlations (hours to seconds) and added memory-safe chord diagram edge capping to prevent OOM crashes |

---

## 6. Limitations & Next Steps

1. **Pathway-Level Biological Annotation:** Map the prioritized KOs (`KO_FC_group.csv`, `differential_functional_contributions.csv`, and `cluego_gene_list_KO.txt`) to specific KEGG/MetaCyc pathways (e.g., GABA/glutamate/tryptophan metabolism, butyrate synthesis, inflammatory signaling) to sharpen mechanistic hypotheses.
2. **Batch & Covariate Sensitivity:** Compare filtered vs. unfiltered `db-RDA` runs and verify consistency across cohorts (`study_id`), age, sex, and hypervariable regions (V3–V4 vs. V4).
3. **Cohort Expansion & Shotgun Integration:** Incorporate additional 16S cohorts using `0_preprocessing/` and integrate clinical subtyping (e.g., ISI severity or insomnia subtypes) and shotgun metagenomic data where available to move from predicted KOs to direct gene/strain resolution.