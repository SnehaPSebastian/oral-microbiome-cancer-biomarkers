# ==========================
# QIIME2 + DADA2 Cheat Sheet
# ==========================

# 1. Data Collection
prefetch SRR9277638 -O sra

# 2. SRA to FASTQ
fasterq-dump sra/SRR9277638/SRR9277638.sra --split-files
for srafile in sra/SRR*/SRR*.sra; do
  fasterq-dump "$srafile" --split-files -O fastq
done

# 3. Check and compress
ls fastq
ls fastq | wc
gzip fastq/*.fastq
ls fastq | grep -v '\.gz$'
rm -i fastq/SRR9277643_*.fastq.gz

# 4. Activate QIIME2
source /opt/homebrew/anaconda3/bin/activate
conda activate qiime2-amplicon-2025.4

# 5. Manifest creation
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
for f in fastq/*_1.fastq.gz; do
  sample=$(basename "$f" _1.fastq.gz)
  echo -e "${sample}\t$(pwd)/fastq/${sample}_1.fastq.gz\t$(pwd)/fastq/${sample}_2.fastq.gz"
done >> manifest.tsv

# 6. Import to QIIME2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 7. Demux Summary
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux-summary.qzv
qiime tools view demux-summary.qzv

# 8. Cutadapt trimming
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-match-adapter-wildcards \
  --p-discard-untrimmed \
  --o-trimmed-sequences demux-trimmed.qza \
  --verbose

# Optional reverse complement
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GGATTAGATACCCBDGTAGTC \
  --p-match-adapter-wildcards \
  --p-discard-untrimmed \
  --o-trimmed-sequences demux-trimmed-test.qza \
  --verbose

qiime demux summarize \
  --i-data demux-trimmed.qza \
  --o-visualization demux-trimmed.qzv

# 9. DADA2 denoise paired
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-trimmed.qza \
  --p-trunc-len-f 270 \
  --p-trunc-len-r 260 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-trimmed.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 225 \
  --p-max-ee-f 3 \
  --p-max-ee-r 4 \
  --p-min-overlap 12 \
  --p-min-fold-parent-over-abundance 4 \
  --p-n-threads 0 \
  --o-table table-v4.qza \
  --o-representative-sequences rep-seqs-v4.qza \
  --o-denoising-stats stats-v4.qza

# Single-end forward
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_fwd_trimmed.qza \
  --p-trunc-len 220 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

# 10. QC & visualization
qiime tools peek table.qza
qiime tools peek rep-seqs.qza
qiime tools peek denoising-stats.qza

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file sample-metadata.tsv \
  --o-visualization table.qzv

# 11. Phylogeny & Diversity
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

# Alpha diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization observed-features-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization evenness-group-significance.qzv

# Beta diversity
qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization bray-curtis-emperor.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column isolate \
  --p-method permanova \
  --p-permutations 999 \
  --o-visualization bray-curtis-permanova.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column isolate \
  --p-method permanova \
  --p-permutations 999 \
  --o-visualization unweighted-unifrac-permanova.qzv

# Export feature table
qiime tools export --input-path table-v4.qza --output-path exported-table-v4
biom convert -i exported-table-v4/feature-table.biom -o feature-table-v4.tsv --to-tsv

# 12. Taxonomy
wget https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza

qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 400 \
  --p-max-length 500 \
  --o-reads ref-seqs-341F-785R.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-341F-785R.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier v3v4-paired-end-classifier.qza

qiime feature-classifier classify-sklearn \
  --i-classifier v3v4-paired-end-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-barplot.qzv

# 13. ANCOM-BC
qiime taxa collapse \
  --i-table table-v4.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table table-l7.qza

qiime composition ancombc \
  --i-table table-l7.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula "isolate" \
  --o-differentials ancombc-results.qza

qiime composition da-barplot \
  --i-data ancombc-results.qza \
  --p-significance-threshold 0.05 \
  --o-visualization ancombc-barplot.qzv

qiime taxa collapse \
  --i-table table-v4.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-l6.qza

qiime composition ancombc \
  --i-table table-l6.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-formula "isolate" \
  --o-differentials l6-ancombc-results.qza

qiime composition da-barplot \
  --i-data l6-ancombc-results.qza \
  --p-significance-threshold 0.05 \
  --p-level-delimiter ';' \
  --o-visualization l6-da-barplot.qzv

qiime composition tabulate \
  --i-data l6-ancombc-results.qza \
  --o-visualization l6-ancombc-table.qzv

