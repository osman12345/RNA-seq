# RNA-Seq Analysis

## Overview

This repository contains an R script for performing RNA-Seq data analysis, including data preprocessing, differential expression analysis, visualization, and functional enrichment. The script utilizes popular R packages such as `edgeR`, `limma`, and `org.Mm.eg.db` for efficient and comprehensive analysis.

---

## Features

### 1. **Data Preprocessing**
- Reads count data and sample information.
- Filters low-expressed genes and normalizes counts.

### 2. **Exploratory Data Analysis**
- Library size bar plots.
- Multidimensional scaling (MDS) plots.
- Heatmaps of highly variable genes.
- Boxplots of log CPM values.

### 3. **Differential Expression Analysis**
- Linear modeling using `voom` from `limma`.
- Generalized linear models (GLM) using `edgeR`.
- Identification of significant differentially expressed genes (DEGs).

### 4. **Visualization**
- MD plots for DEGs.
- Heatmaps of top DEGs.
- Venn diagrams for comparison of conditions.

### 5. **Functional Enrichment**
- Gene Ontology (GO) analysis.
- KEGG pathway analysis.
- Top enriched pathways and biological processes.

---

## Dependencies

The following R packages are required to run the script:

- [`edgeR`](https://bioconductor.org/packages/edgeR)
- [`limma`](https://bioconductor.org/packages/limma)
- [`org.Mm.eg.db`](https://bioconductor.org/packages/org.Mm.eg.db)
- [`RColorBrewer`](https://cran.r-project.org/package=RColorBrewer)
- [`gplots`](https://cran.r-project.org/package=gplots)
- [`tidyverse`](https://cran.r-project.org/package=tidyverse)

To install these packages, run:
```R
install.packages(c("RColorBrewer", "gplots", "tidyverse"))
BiocManager::install(c("edgeR", "limma", "org.Mm.eg.db"))
