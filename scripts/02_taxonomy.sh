#!/usr/bin/env bash
# QIIME2 Taxonomic Classification
# Step 02: Assign taxonomy and generate taxa bar plots
# Usage: bash scripts/02_taxonomy.sh

set -euo pipefail

# -------------------------------
# 1. Define input/output paths
# -------------------------------
TABLE="results/table.qza"
REP_SEQS="results/rep-seqs.qza"
METADATA="data/metadata.tsv"

# Reference classifier
CLASSIFIER="reference/silva-138-99-nb-classifier.qza"

# Output folder
TAX_DIR="results/taxonomy"
mkdir -p "$TAX_DIR"

TAXONOMY="$TAX_DIR/taxonomy.qza"
TAXONOMY_VIZ="$TAX_DIR/taxonomy.qzv"
TAXA_BAR="$TAX_DIR/taxa-bar-plots.qzv"

# -------------------------------
# 2. Taxonomic classification
# -------------------------------
echo ">>> Assigning taxonomy..."
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$REP_SEQS" \
  --o-classification "$TAXONOMY"

# -------------------------------
# 3. Visualize taxonomy
# -------------------------------
echo ">>> Creating taxonomy visualization..."
qiime metadata tabulate \
  --m-input-file "$TAXONOMY" \
  --o-visualization "$TAXONOMY_VIZ"

# -------------------------------
# 4. Taxa bar plots
# -------------------------------
echo ">>> Generating taxa bar plots..."
qiime taxa barplot \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$TAXA_BAR"

echo ">>> Step 02 complete: Taxonomy assigned!"

