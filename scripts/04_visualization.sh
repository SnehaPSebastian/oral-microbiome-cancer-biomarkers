#!/bin/bash
# Visualization and export workflow

# Step 0: Validate inputs
echo "Validating taxonomy..."
qiime tools validate results/taxonomy.qza

echo "Validating feature table..."
qiime tools validate results/table.qza

# Step 1: Taxa bar plots (already created, but rerun if needed)
qiime taxa barplot \
  --i-table results/table.qza \
  --i-taxonomy results/taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization results/taxa-bar-plots.qzv

# Step 2: Collapse table at genus level
qiime taxa collapse \
  --i-table results/table.qza \
  --i-taxonomy results/taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table results/genus-table.qza

# Step 3: Export collapsed table to TSV
qiime tools export \
  --input-path results/genus-table.qza \
  --output-path exports/genus-table

# Step 4: Convert BIOM to TSV
biom convert \
  -i exports/genus-table/feature-table.biom \
  -o exports/genus-table.tsv \
  --to-tsv
