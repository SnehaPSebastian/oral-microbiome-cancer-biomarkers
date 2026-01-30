#Statistical analysis of differential taxa (R). v4

#Setting directory
setwd("C:/Users/Sneha P Sebastian/OneDrive/Desktop/Fastq/Oral_microbiome")

#Installing and loading the packages
#install.packages(c("tidyverse", "reshape2"))
library(tidyverse)
library(reshape2)

#Loading feature table 
otu <- read.delim(
  "feature-table-v4.tsv",
  row.names = 1,
  skip = 1,
  check.names = FALSE
)
otu #To view the table

#Confirmation checks
colnames(otu)[1:5]#valid sample IDs.
dim(otu)          # how big the table is
summary(otu[,1]) # check one sample

#Loading metadata
meta <- read.delim(
  "Metadata_toba.tsv",
  row.names = 1,
  check.names = FALSE
)
meta

#Remove low depth samples
sample_sums <- colSums(otu)
otu <- otu[, sample_sums > 1000]
summary(sample_sums)


#match samples
meta <- meta[colnames(otu), ]
any(is.na(meta)) #checking

#Normalising the data
otu_rel <- sweep(otu, 2, colSums(otu), "/")
colSums(otu_rel)[1:5] #Cheking if normalisation worked

#Preparing grouping variables
colnames(meta) #Checking columns
table(meta$Group)

#Removing rare taxa
prev <- rowSums(otu_rel > 0)
otu_filt <- otu_rel[prev >= 0.1 * ncol(otu_rel), ] #This removes noisy, extremely rare ASVs

# convert the meta$Group to a factor
meta$Group <- factor(meta$Group,
                     levels = c("Control", "TobaccoAbuser", "OralCancer"))
# Check
table(meta$Group)

#Kruskal Wallis test for one taxon (example)
taxon1 <- otu_rel[1, ]

df <- data.frame(
  value = as.numeric(taxon1),
  Group = meta$Group
)

kruskal.test(value ~ Group, data = df)


#Kruskal Wallis test for ALL taxa (final version)
pvals <- apply(otu_rel, 1, function(x) {
  df <- data.frame(
    value = as.numeric(x),
    Group = meta$Group
  )
  kruskal.test(value ~ Group, data = df)$p.value
})

pvals_adj <- p.adjust(pvals, method = "fdr")

#Results table 
Kruskal_results <- data.frame(
  Taxon = rownames(otu_rel),
  p_value = pvals,
  p_adj = pvals_adj
)

head(Kruskal_results)

#Find significant taxa
sig_taxa <- Kruskal_results[Kruskal_results$p_adj < 0.05, ]
nrow(sig_taxa)

#Simple visualization
top_taxa <- Kruskal_results$Taxon[1:3]

par(mfrow = c(1,3))
for (tx in top_taxa) {
  df_plot <- data.frame(
    value = as.numeric(otu_filt[tx, ]),
    Group = meta$Group
  )
  
  boxplot(
    value ~ Group,
    data = df_plot,
    main = tx,
    ylab = "Relative abundance",
    xlab = "Group",
    title = "Kruskal Wallis test"
  )
}
par(mfrow = c(1,1))

#########################################################################
#Alpha Diversity Analysis

#Loading libraries
#install.packages("vegan")
library(vegan)
library(ggplot2)

# Transpose: samples as rows
otu_t <- t(otu)

# Shannon diversity
shannon <- diversity(otu_t, index = "shannon")

# Observed richness (number of ASVs)
observed <- specnumber(otu_t)

#Create aplha dataframe
alpha_df <- data.frame(
  SampleID = rownames(otu_t),
  Shannon = shannon,
  Observed = observed,
  Group = meta$Group
)
# Sanity check
head(alpha_df)
table(alpha_df$Group)

#Box plots for Alpha diversity
#Shannon diversity
boxplot(
  Shannon ~ Group,
  data = alpha_df,
  main = "Alpha diversity (Shannon)",
  ylab = "Shannon index",
  xlab = "Group")

#Observed ASVs
boxplot(
  Observed ~ Group,
  data = alpha_df,
  main = "Alpha diversity (Observed ASVs)",
  ylab = "Observed features",
  xlab = "Group")

#Statistical test (Shannon and Observed)
kruskal.test(Shannon ~ Group, data = alpha_df)
kruskal.test(Observed ~ Group, data = alpha_df)

#Pairwise test (Since Shannon is significant)
pairwise.wilcox.test(
  alpha_df$Shannon,
  alpha_df$Group,
  p.adjust.method = "fdr")

##########################################################################
# Beta diversity Analysis in R

#Import Libraries
library(vegan)
library(ggplot2)

#Preparing OTU table with relative abundance data
# Transpose: samples as rows
otu_rel_t <- t(otu_rel)
# Sanity check
dim(otu_rel_t)

#Calculating the Bray-Curtis Distance
bray_dist <- vegdist(otu_rel_t, method = "bray")
# Check
bray_dist

#PCoA (Principal Coordinate Analysis)
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
# Create dataframe
pcoa_df <- data.frame(
  SampleID = rownames(pcoa$points),
  PC1 = pcoa$points[,1],
  PC2 = pcoa$points[,2],
  Group = meta$Group)
# % variance explained
var_exp <- round(100 * pcoa$eig / sum(pcoa$eig), 1)

#PCoA plot
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCoA – Bray-Curtis distance",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 12))

#Pairwise PERMANOVA
adonis2(bray_dist ~ Group, data = meta, permutations = 999)

groups <- levels(meta$Group)
pairwise_results <- data.frame(
  Comparison = character(),
  F_model = numeric(),
  R2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    
    keep <- meta$Group %in% c(groups[i], groups[j])
    
    dist_sub <- as.dist(as.matrix(bray_dist)[keep, keep])
    meta_sub <- meta[keep, , drop = FALSE]
    
    res <- adonis2(dist_sub ~ Group, data = meta_sub, permutations = 999)
    
    pairwise_results <- rbind(
      pairwise_results,
      data.frame(
        Comparison = paste(groups[i], "vs", groups[j]),
        F_model = res$F[1],
        R2 = res$R2[1],
        p_value = res$`Pr(>F)`[1]
      )
    )
  }
}

# FDR correction
pairwise_results$p_adj <- p.adjust(pairwise_results$p_value, method = "fdr")
pairwise_results

#Checking dispersion to validate PERMANOVA
disp <- betadisper(bray_dist, meta$Group)
anova(disp)


#Rarefaction depth

set.seed(123)

otu_rare <- rrarefy(otu_t, sample = min_depth) #Randomly subsamples exactly 10,010 reads from each sample
                                               #Without replacement

dim(otu_rare) #TO check the dimetions
rowSums(otu_rare) #computes total reads in each sample:

min_depth <- min(rowSums(otu_t)) #Find the smallest library size among all samples
cat("The rarefraction depth:\n",min_depth)


shannon_rare <- diversity(otu_rare, index = "shannon")
observed_rare <- specnumber(otu_rare)

kruskal.test(shannon_rare ~ meta$Group)
kruskal.test(observed_rare ~ meta$Group)


#Relative abundance

otu_rel <- sweep(otu_t, 1, rowSums(otu_t), FUN = "/") #converts each sample to proportions
                                                      #each row now sums to 1
                                                      #Values represent fraction of community
rowSums(otu_rel) #to verify


otu_rel[1, 1] * 100 #relative abundance for one sample and converting it to percentage

#Add Group info to OTU table
otu_rel_df <- as.data.frame(otu_rel)
otu_rel_df$Group <- meta$Group

#Calculate mean relative abundance per ASV per group
library(dplyr)

group_means <- otu_rel_df %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), mean))

#Extract Top 10 ASVs per group
top10_long <- group_means %>%
  pivot_longer(
    cols = -Group,
    names_to = "ASV",
    values_to = "Mean_RelAbundance"
  ) %>%
  group_by(Group) %>%
  arrange(desc(Mean_RelAbundance)) %>%
  slice_head(n = 10)
top10_long

#Attaching taxonomy to ASVs
#Load the taxonomy
tax <- read.table(
  "taxonomy.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
head(tax)

#Split taxonomy into levels & clean prefixes
library(dplyr)
library(tidyr)

tax_clean <- tax %>%     #removes d__  p__  g__
  separate(
    Taxon,
    into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
    sep = "; ",
    fill = "right"
  ) %>%
  mutate(across(
    Domain:Species,
    ~ gsub("^[a-z]__", "", .)
  ))
head(tax_clean)

#Merge taxonomy with relative abundance
otu_long <- otu_rel_df %>%      #Converting to long format
  pivot_longer(
    cols = -Group,
    names_to = "Feature.ID",
    values_to = "RelAbundance"
  )

#Merging
otu_tax <- otu_long %>%
  left_join(tax_clean, by = "Feature.ID")

#Collapse to GENUS level (key step)
genus_rel <- otu_tax %>%
  group_by(Group, Genus) %>%
  summarise(
    Mean_RelAbundance = mean(RelAbundance),
    .groups = "drop"
  )
head(genus_rel)

#Extract Top 10 genera per group
top10_genus <- genus_rel %>%
  group_by(Group) %>%
  arrange(desc(Mean_RelAbundance)) %>%
  slice_head(n = 10)
top10_genus

#Stacked bar plot (Top genera by group)
library(ggplot2)

top10_genus
#Convert to percentages
top10_genus <- top10_genus %>%
  mutate(Percent = Mean_RelAbundance * 100)

#Plot the stacked bar plot
ggplot(top10_genus, aes(x = Group, y = Percent, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Top 10 Genera by Relative Abundance",
    y = "Mean Relative Abundance (%)",
    x = "Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Phylum level

# ============================================================
# Microbiome Relative Abundance Analysis (OTU/Taxonomy/Metadata)
# ============================================================

# 1️⃣ Load libraries
library(qiime2R)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ========================
# 2️⃣ Paths to your files
# ========================
otu_path <- "C:/Users/Sneha P Sebastian/OneDrive/Desktop/Fastq/Oral_microbiome/table-v4.qza"
tax_path <- "C:/Users/Sneha P Sebastian/OneDrive/Desktop/Fastq/Oral_microbiome/taxonomy.qza"
meta_path <- "C:/Users/Sneha P Sebastian/OneDrive/Desktop/Fastq/Oral_microbiome/Metadata_toba.tsv"

# ========================
# 3️⃣ Import OTU table
# ========================
otu_qza <- read_qza(otu_path)
otu_table <- as.data.frame(otu_qza$data)

# Transpose if OTUs are columns
if(nrow(otu_table) < ncol(otu_table)){
  otu_table <- t(otu_table)
}

# Add OTU IDs if missing
if(is.null(rownames(otu_table))){
  rownames(otu_table) <- paste0("OTU", seq_len(nrow(otu_table)))
}

# ========================
# 4️⃣ Import taxonomy table
# ========================
tax_qza <- read_qza(tax_path)
tax_table <- as.data.frame(tax_qza$data)

# Split Taxon column into ranks if present
if("Taxon" %in% colnames(tax_table)){
  tax_split <- as.data.frame(do.call(rbind, strsplit(as.character(tax_table$Taxon), "; ")))
  colnames(tax_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_split <- as.data.frame(lapply(tax_split, function(x) gsub("^[a-z]__","", x)))
  tax_split[is.na(tax_split)] <- "Unknown"
  tax_split[tax_split==""] <- "Unknown"
  rownames(tax_split) <- rownames(otu_table)
  tax_table <- tax_split
}

# ========================
# 5️⃣ Import metadata
# ========================
sample_data <- read.delim(meta_path, row.names = 1)

# Keep only shared samples
shared_samples <- intersect(colnames(otu_table), rownames(sample_data))
otu_table <- otu_table[, shared_samples]
sample_data <- sample_data[shared_samples, ]

# ========================
# 6️⃣ Create phyloseq object
# ========================
ps <- phyloseq(
  otu_table(otu_table, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_table)),
  sample_data(sample_data)
)

ps  # Check the object

# ========================
# 7️⃣ Transform to relative abundance
# ========================
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# ========================
# 8️⃣ Aggregate by taxonomic level
# ========================
ps_phylum  <- tax_glom(ps_rel, taxrank = "Phylum")
ps_genus   <- tax_glom(ps_rel, taxrank = "Genus")
ps_species <- tax_glom(ps_rel, taxrank = "Species")

# ========================
# 9️⃣ Melt for ggplot2
# ========================
phylum_df  <- psmelt(ps_phylum)
genus_df   <- psmelt(ps_genus)
species_df <- psmelt(ps_species)

# ========================
# 🔟 Keep top N taxa, collapse others
# ========================
top_n <- 10

collapse_low <- function(df, tax_col){
  top_taxa <- df %>%
    group_by_at(tax_col) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
    top_n(top_n, mean_abundance) %>%
    pull(!!sym(tax_col))
  
  df[[tax_col]] <- as.character(df[[tax_col]])
  df[[tax_col]][!(df[[tax_col]] %in% top_taxa)] <- "Other"
  return(df)
}

phylum_df  <- collapse_low(phylum_df, "Phylum")
genus_df   <- collapse_low(genus_df, "Genus")
species_df <- collapse_low(species_df, "Species")

# ========================
# 1️⃣1️⃣ Plot function
# ========================
# Top 10 Phyla by mean relative abundance
phylum_top10 <- phylum_df %>%
  group_by(Group, Phylum) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  ungroup() %>%
  group_by(Phylum) %>%
  summarise(overall_mean = mean(mean_abundance)) %>%
  slice_max(order_by = overall_mean, n = 10) %>%
  pull(Phylum)

phylum_plot_df <- phylum_df %>%
  filter(Phylum %in% phylum_top10) %>%
  group_by(Group, Phylum) %>%
  summarise(mean_abundance = mean(Abundance))

ggplot(phylum_plot_df,
       aes(x = Group, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab("Mean Relative Abundance (%)") +
  xlab("Group") +
  ggtitle("Top 10 Phyla by Relative Abundance")

# Remove unclassified species (recommended)
species_df_clean <- species_df %>%
  filter(!is.na(Species), Species != "NA")

# Top 10 Species by mean relative abundance
species_df <- psmelt(ps_species) %>%
  dplyr::select(Sample, Group, Species, Abundance)
head(species_df)

species_df_clean <- species_df %>%
  filter(!is.na(Species),
         Species != "",
         !grepl("uncultured|unidentified", Species, ignore.case = TRUE))

species_top10 <- species_df_clean %>%
  group_by(Group, Species) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  group_by(Species) %>%
  summarise(overall_mean = mean(mean_abundance)) %>%
  slice_max(overall_mean, n = 10) %>%
  pull(Species)

species_plot_df <- species_df_clean %>%
  filter(Species %in% species_top10) %>%
  group_by(Group, Species) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")
ggplot(species_plot_df,
       aes(x = Group, y = mean_abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ylab("Mean Relative Abundance (%)") +
  xlab("Group") +
  ggtitle("Top 10 Species by Relative Abundance")


