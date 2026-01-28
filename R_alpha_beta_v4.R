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

#Rarefraction

set.seed(123)
min_depth <- min(rowSums(otu_t))
otu_rare <- rrarefy(otu_t, sample = min_depth)
dim(otu_rare)
rowSums(otu_rare)

shannon_rare <- diversity(otu_rare, index = "shannon")
observed_rare <- specnumber(otu_rare)

kruskal.test(shannon_rare ~ meta$Group)
kruskal.test(observed_rare ~ meta$Group)

