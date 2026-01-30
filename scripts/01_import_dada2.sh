#!/usr/bin/env bash
# QIIME2 Import + DADA2 Pipeline
# Step 01: Import FASTQ files and run DADA2 denoising
# Usage: bash scripts/01_import_dada2.sh

set -euo pipefail   # safer bash: stop on errors, undefined vars

# -------------------------------
# 1. Define input/output paths
# -------------------------------
MANIFEST="data/manifest.tsv"
METADATA="data/metadata.tsv"
DEMUX="results/paired-end-demux.qza"
DEMUX_SUMMARY="results/demux-summary.qzv"
TRIMMED="results/demux-trimmed.qza"
TABLE="results/table.qza"
REP_SEQS="results/rep-seqs.qza"
STATS="results/denoising-stats.qza"

# -------------------------------
# 2. Import FASTQ files
# -------------------------------
echo ">>> Importing FASTQ files using manifest..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path "$DEMUX" \
  --input-format PairedEndFastqManifestPhred33V2

# -------------------------------
# 3. Summarize demultiplexed data
# -------------------------------
echo ">>> Summarizing demultiplexed sequences..."
qiime demux summarize \
  --i-data "$DEMUX" \
  --o-visualization "$DEMUX_SUMMARY"

# -------------------------------
# 4. Adapter trimming (Cutadapt)
# -------------------------------
echo ">>> Trimming adapters..."
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "$DEMUX" \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-match-adapter-wildcards \
  --p-discard-untrimmed \
  --o-trimmed-sequences "$TRIMMED" \
  --verbose

# -------------------------------
# 5. DADA2 denoising
# -------------------------------
echo ">>> Running DADA2 denoise-paired..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$TRIMMED" \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 225 \
  --o-table "$TABLE" \
  --o-representative-sequences "$REP_SEQS" \
  --o-denoising-stats "$STATS"

# -------------------------------
# 6. Quick checks
# -------------------------------
echo ">>> Checking outputs..."
qiime tools peek "$TABLE"
qiime tools peek "$REP_SEQS"
qiime tools peek "$STATS"

echo ">>> Step 01 complete: Import + DADA2 finished!"
