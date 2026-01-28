# Exploring the Oral Microbiome for Biomarker Discovery in Early Cancer Detection

## Project Overview
This repository documents the use of Linux command-line tools and Git/GitHub version control
for managing a bioinformatics project focused on the oral microbiome.  
The aim is to understand basic workflow organization, reproducibility, and collaboration
using Git while working with microbiome-related data and analysis scripts.

## Objectives
- Learn and practice Linux command-line operations
- Understand Git version control and GitHub repositories
- Organize a bioinformatics project directory
- Maintain reproducible analysis using scripts and documentation
- Apply these skills to an oral microbiome research context

## Directory Structure
Oral_microbiome/
├── scripts/ # Analysis scripts
├── Exploring-the-Oral-Microbiome-for-Biomarker-Discovery-in-Early-Cancer-Detection/
│ └── documentation and notes
├── R_alpha_beta_v4.R # R script for diversity analysis
├── .gitignore # Files excluded from version control
└── README.md


## Tools and Environment
- Operating System: Linux (WSL)
- Version Control: Git and GitHub
- Programming Language: R
- Bioinformatics Framework: QIIME2 (used locally)

## Git Workflow Used
The following Git commands were used during this project:

```bash
git init
git branch -m main
git add <files>
git commit -m "Initial commit"
git remote add origin <repository_url>
git push -u origin main

After setting the upstream branch, subsequent updates were pushed using:

git push

Data Management

Raw sequencing data and large intermediate files were not uploaded to GitHub due to size
constraints and best practices.

The following files and directories are excluded using .gitignore:

Raw FASTQ files

QIIME2 artifacts (.qza, .qzv)

Large processed tables (.tsv, .biom)

System and environment files

All analysis steps are documented through scripts and notes to ensure reproducibility.

Collaboration

This repository is designed for collaborative work using GitHub.
Contributors can clone the repository, create branches, and push changes following standard Git workflows.

Author(s)

Aisha Hassan Blahayil
Sneha P Sebastian


