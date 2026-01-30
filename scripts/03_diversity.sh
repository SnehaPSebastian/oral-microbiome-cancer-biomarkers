#!/bin/bash
# Diversity analysis workflow

# Step 0: Validate inputs
echo "Validating feature table..."
qiime tools validate results/table.qza

echo "Validating taxonomy..."
qiime tools validate results/taxonomy.qza

echo "Validating metadata..."
test -f metadata.tsv || { echo "metadata.tsv not found!"; exit 1; }

# Step 1: Core diversity metrics (requires a phylogenetic tree if you want Faith’s PD)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny results/rooted-tree.qza \
  --i-table results/table.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file metadata.tsv \
  --output-dir results/core-metrics-results

# Step 2: Alpha diversity visualizations
qiime diversity alpha-group-significance \
  --i-alpha-diversity results/core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization results/core-metrics-results/shannon-group-significance.qzv

# Step 3: Beta diversity visualizations
qiime diversity beta-group-significance \
  --i-distance-matrix results/core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column SampleType \
  --o-visualization results/core-metrics-results/bray-curtis-group-significance.qzv \
  --p-pairwise

# Step 4: Ordination (PCoA + Emperor)
qiime emperor plot \
  --i-pcoa results/core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization results/core-metrics-results/bray-curtis-emperor.qzv
