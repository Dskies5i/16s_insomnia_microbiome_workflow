#!/usr/bin/env bash
#NAME: 1_meta_analysis/10A_run_picrust2.sh
#COMMENT: PICRUSt2 runner using SEPP placement for macOS/Linux stability.
#USE: conda activate qiime2-amplicon-2024.5
#     bash 1_meta_analysis/10A_run_picrust2.sh

set -euo pipefail

echo "===================================================="
echo "  PICRUSt2 AUTOMATED EXECUTION PIPELINE (Script 10A)"
echo "===================================================="
echo ""

# 1. Safety check: must be in the correct conda environment
if [[ "${CONDA_DEFAULT_ENV:-}" != "qiime2-amplicon-2024.5" ]]; then
    echo "❌ Error: This script must be run from the 'qiime2-amplicon-2024.5' conda environment."
    echo "Run: conda activate qiime2-amplicon-2024.5"
    exit 1
fi
echo "✔ Conda environment verified: $(conda info --envs | grep '*' | awk '{print $1}')"

# 2. Install PICRUSt2 if not present in the environment
if ! command -v picrust2_pipeline.py &>/dev/null; then
    echo "⚙ PICRUSt2 not found in current environment. Installing via Conda..."
    conda install -c bioconda -c conda-forge picrust2 -y
else
    echo "✔ PICRUSt2 is already installed."
fi

# 3. Dynamic path resolution from config.yaml
echo "🔍 Extracting input paths from config/config.yaml..."
if [ ! -f "config/config.yaml" ]; then
    echo "❌ Error: config/config.yaml not found in root directory."
    exit 1
fi

TSV=$(Rscript --vanilla -e "suppressWarnings(suppressPackageStartupMessages(cfg <- yaml::read_yaml('config/config.yaml'))); cat(cfg\$input\$otu)" 2>/dev/null | tr -d '\r')

if [ -z "$TSV" ]; then
    echo "⚠ Warning: Dynamic config extraction failed. Falling back to default path."
    TSV="data/pre_processed/unified/unified_otu_table.tsv"
fi

DATA_DIR=$(dirname "$TSV")
FASTA="$DATA_DIR/unified_rep_seqs.fasta"
BIOM="$DATA_DIR/feature-table.biom"

echo "  - Resolved OTU Table (TSV): $TSV"
echo "  - Resolved Representative Sequences (FASTA): $FASTA"
echo "  - Target BIOM path: $BIOM"

# 4. Check and prepare required input files
if [ ! -f "$FASTA" ]; then
    echo "❌ Error: Missing required input FASTA file at: $FASTA"
    exit 1
fi

if [ ! -f "$BIOM" ]; then
    if [ -f "$TSV" ]; then
        echo "⚠ feature-table.biom not found, but unified_otu_table.tsv exists."
        echo "⚙ Converting TSV to BIOM format..."
        biom convert -i "$TSV" -o "$BIOM" --table-type="OTU table"
        echo "✔ Successfully converted TSV to BIOM."
    else
        echo "❌ Error: Neither feature-table.biom nor unified_otu_table.tsv found."
        exit 1
    fi
fi

echo "✔ Input files verified successfully."
echo ""

# 5. Sanity check: BIOM should contain samples
NUM_SAMPLES=$(biom summarize-table -i "$BIOM" | grep "Num samples:" | awk '{print $3}' | tr -d ',')
echo "📊 BIOM sample count: $NUM_SAMPLES"
if [ "$NUM_SAMPLES" -le 1 ]; then
    echo "❌ Error: BIOM reports 0 or 1 samples. Refusing to run."
    exit 1
fi

# 6. Checkpoint: Skip computation if output already exists
FINAL_OUTPUT="outputs/tables/merged_KO_metagenome.tsv"
if [ -f "$FINAL_OUTPUT" ]; then
    echo ""
    echo "⏭ [SALTADO] PICRUSt2 output already exists at: $FINAL_OUTPUT"
    echo " (Se omite la computación pesada y se continúa con el pipeline)."
    echo ""
    exit 0
fi

# 7. Run PICRUSt2 Pipeline with SEPP Placement (Bypasses epa-ng C++ crashes)
echo ""
echo "🚀 Running PICRUSt2 pipeline with SEPP placement..."
echo ""

if [ -d "picrust2_out" ]; then
    rm -rf picrust2_out
fi

ulimit -s unlimited 2>/dev/null || true

if [[ "$(uname)" == "Darwin" ]]; then
    THREADS=6
    echo "🍏 Detected macOS: setting threads to $THREADS."
else
    THREADS=10
    echo "🐧 Detected Linux: setting threads to $THREADS."
fi

set +e
picrust2_pipeline.py \
    -s "$FASTA" \
    -i "$BIOM" \
    -o picrust2_out \
    --placement_tool sepp \
    -p "$THREADS"
PIPELINE_EXIT_CODE=$?
set -e

if [ $PIPELINE_EXIT_CODE -ne 0 ]; then
    echo "❌ Error: PICRUSt2 pipeline failed with exit code $PIPELINE_EXIT_CODE."
    exit $PIPELINE_EXIT_CODE
fi

# 8. Locate, extract, and format the output KO table
echo ""
echo "📦 Formatting and moving output files..."
mkdir -p outputs/tables

FOUND_KO=$(find picrust2_out -name "pred_metagenome_unstrat.tsv.gz" | head -n 1)

if [ -z "$FOUND_KO" ] || [ ! -f "$FOUND_KO" ]; then
    echo "❌ Error: Could not find pred_metagenome_unstrat.tsv.gz inside picrust2_out/"
    exit 1
fi

echo "✔ Found KO prediction file at: $FOUND_KO"
cp "$FOUND_KO" outputs/tables/pred_metagenome_unstrat.tsv.gz
gunzip -f outputs/tables/pred_metagenome_unstrat.tsv.gz
mv outputs/tables/pred_metagenome_unstrat.tsv "$FINAL_OUTPUT"

echo ""
echo "✔ Successfully created $FINAL_OUTPUT!"
echo ""
echo "First 2 rows of the generated KO table:"
head -2 "$FINAL_OUTPUT"
echo ""

# FUNCIONA PERO ABORTA LA PIPELINE... 
## 9. Verify Sample ID Alignment against Phyloseq
# echo "===================================================="
# echo "  VERIFYING SAMPLE ID ALIGNMENT"
# echo "===================================================="
# echo "First 3 Sample IDs in merged_KO_metagenome.tsv:"
# head -1 "$FINAL_OUTPUT" | tr '\t' '\n' | sed -n '2,4p'
# 
# if [ -f "outputs/models/physeq_filtered.rds" ]; then
#     echo ""
#     echo "First 3 Sample IDs in physeq_filtered.rds:"
#     Rscript -e 'ps <- readRDS("outputs/models/physeq_filtered.rds"); cat(head(sample_names(ps), 3), sep = "\n")'
# else
#     echo "⚠ Warning: outputs/models/physeq_filtered.rds not found yet to cross-check sample IDs."
# fi

# 9. Verify Sample ID Alignment against Phyloseq
echo "==================================================--"
echo "  VERIFYING SAMPLE ID ALIGNMENT"
echo "==================================================--"
echo "First 3 Sample IDs in merged_KO_metagenome.tsv:"
head -1 "$FINAL_OUTPUT" | tr '\t' '\n' | sed -n '2,4p'

if [ -f "outputs/models/physeq_filtered.rds" ]; then
    echo ""
    echo "First 3 Sample IDs in physeq_filtered.rds:"
    Rscript -e 'suppressPackageStartupMessages(library(phyloseq)); ps <- readRDS("outputs/models/physeq_filtered.rds"); cat(head(sample_names(ps), 3), sep = "\n")'
else
    echo "⚠ Warning: outputs/models/physeq_filtered.rds not found yet to cross-check sample IDs."
fi

echo ""
echo "🎉 PICRUSt2 execution complete! You can now proceed to script 10B."
echo "Rscript 1_meta_analysis/10B_function_taxa_ko_cor.R"
echo ""

echo ""
echo "🎉 PICRUSt2 execution complete! You can now proceed to script 10B."
echo "Rscript 1_meta_analysis/10B_function_taxa_ko_cor.R"
echo ""




# 
# 
# 
# 
# 
#
#
#
#
