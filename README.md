# Exploring the Oral Microbiome for Biomarker Discovery in Early Cancer Detection

## Overview

This repository documents the **Linux, Git, and microbiome bioinformatics skills** acquired as part of Mini project. The work focuses on understanding and implementing a **QIIME 2–style amplicon sequencing analysis workflow**, using modern QIIME2 tools while retaining the original conceptual logic of QIIME 1.

**Important note**: This repository contains **only scripts, commands, and documentation**. Large biological datasets and QIIME-generated artifacts are intentionally excluded via `.gitignore` to follow best practices in reproducible computational biology.

The emphasis is on workflow design, command-line proficiency, and version control rather than raw data storage.

---

## Learning Objectives

* Develop proficiency in **Linux command-line environments** for biological data analysis
* Apply **Git version control** for tracking, collaboration, and reproducibility
* Understand the end-to-end **QIIME microbiome analysis workflow** for 16S rRNA data
* Implement **best practices in computational reproducibility and project organization****

---

## Linux Basics Used

The following Linux commands were practiced during the assignment:

```bash
pwd                # check current directory
ls                 # list files
cd directory_name  # navigate directories
mkdir data scripts # create directories
touch file.txt     # create files
rm file.txt        # remove files
cp source dest     # copy files
mv old new         # rename/move files
```

These commands were used to organize sequencing data, scripts, and results.

---

## Git & GitHub Workflow

### Initial Setup

```bash
git init
git branch -m main
git remote add origin <repository-url>
```

### Daily Workflow

```bash
git status
git add filename
git commit -m "Meaningful message"
git push
```

### Best Practices Followed

* Large files excluded using `.gitignore`
* Scripts and documentation tracked
* Repository kept lightweight and reproducible

---

## QIIME Workflow (QIIME 1 Logic)

Although QIIME2 was used, the **conceptual workflow follows QIIME 1 principles**.

### 1. Input Data

* Paired-end FASTQ files
* Sample metadata (`sample-metadata.tsv`)

*(Raw FASTQ files are not included in this repository)*

---

### 2. Demultiplexing & Quality Control

Purpose: Assess read quality and remove low-quality sequences.

Typical steps:

* Demultiplex sequences
* Trim adapters and primers
* Visualize quality profiles

Outputs:

* Demultiplexed reads
* Quality summary visualizations

---

### 3. Feature Table Construction

Purpose: Convert reads into representative sequence features.

Steps:

* Denoising / OTU-style clustering
* Generate feature table
* Generate representative sequences

Outputs:

* Feature table
* Representative sequences

---

### 4. Phylogenetic Tree Construction

Purpose: Enable phylogenetic diversity metrics.

Steps:

* Sequence alignment
* Masking hypervariable regions
* Tree construction (rooted and unrooted)

---

### 5. Diversity Analysis

#### Alpha Diversity

* Observed features
* Shannon diversity
* Faith’s phylogenetic diversity

#### Beta Diversity

* Bray–Curtis distance
* Unweighted UniFrac
* PERMANOVA for group comparison

These analyses help identify **community-level differences between sample groups**.

---

### 6. Statistical Testing

Purpose: Identify statistically significant differences between groups.

Examples:

* Alpha diversity group significance
* Beta diversity PERMANOVA

---

## What Is NOT Included (By Design)

* Raw sequencing data (FASTQ)
* QIIME artifacts (`.qza`, `.qzv`)
* Feature tables and BIOM files

These files are excluded to:

* Keep the repository lightweight
* Follow GitHub best practices
* Ensure reproducibility via scripts

---

## Repository Structure

```text
oral-microbiome-cancer-biomarkers/
├── README.md
├── scripts/          # analysis scripts
├── .gitignore        # excludes large data files
```

---

## Key Takeaways

* Linux enables efficient data handling
* Git ensures reproducibility and collaboration
* Microbiome workflows follow logical, modular steps
* Good bioinformatics practice separates **code** from **data**

---

## Contributors

* **Sneha P Sebastian**
* **Aisha Hassan Blahayil**

---

## Assignment Reflection

This assignment strengthened practical understanding of how microbiome analyses are structured computationally, from raw sequencing reads to ecological interpretation. Emphasis was placed on separating data from code, documenting workflows clearly, and using version control effectively—skills that are essential for scalable and collaborative bioinformatics research.

---

## Disclaimer

This repository is intended strictly for **academic and educational purposes**. No clinical, patient-identifiable, or sensitive biological data are included.

