# scTidyGene: An R Package for scRNA-seq Gene Expression Analysis

## Description
The `scTidyGene` R package facilitates the automated identification and generation of lists containing Conserved Expressed Genes (CEGs) and Differentially Expressed Genes (DEGs) during single-cell RNA sequencing (scRNA-seq) analysis. The package processes all cell clusters (cell types) in a loop, automatically calculating and generating lists of CEGs and/or DEGs. It includes functionality to detect and exclude clusters with fewer than three cells, providing automated reporting to ensure robust analysis.

## Installation
To install `scTidyGene` from GitHub, use the following commands in R:

```R
install.packages("devtools")  # If not already installed
library(devtools)
devtools::install_github("LingzhangMeng/scTidyGene")
```

## Dependencies
The `scTidyGene` package requires the following R packages:

```R
library(Seurat)
library(progress)
library(dplyr)
library(tidyverse)
library(data.table)
```

Ensure these dependencies are installed prior to using `scTidyGene`.

## User Guide

### Function 1: `scCEGs()`
This function calculates and generates a list of Conserved Expressed Genes (CEGs) across specified groups.

**Usage:**
```R
scCEGs(seu_obj, grouping.var)
```

**Parameters:**
- `seu_obj`: A Seurat object containing at least two integrated groups.
- `grouping.var`: A character string specifying the grouping factor (e.g., "condition").

**Output:**
A list of conserved expressed genes.

### Function 2: `scDEGs()`
This function calculates and generates a list of Differentially Expressed Genes (DEGs) between two specified groups.

**Usage:**
```R
scDEGs(seu_obj, ident.1, ident.2, group.by, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
```

**Parameters:**
- `seu_obj`: A Seurat object containing at least two integrated groups.
- `ident.1`: A character string specifying the first group (numerator) for differential expression analysis.
- `ident.2`: A character string specifying the second group (denominator) for differential expression analysis.
- `group.by`: A character string specifying the grouping factor (e.g., "condition").
- `min.pct`: A numeric value indicating the minimum percentage of cells expressing a gene (default: 0.25).
- `logfc.threshold`: A numeric value specifying the minimum log2 fold change threshold for gene expression (default: 0.25).
- `only.pos`: A logical value. If `TRUE`, returns only DEGs with a p-value < 0.05 (default). If `FALSE`, returns all DEGs, including those with a p-value ≥ 0.05.

**Output:**
A list of differentially expressed genes based on the specified parameters.

## Example Workflow

### Data Preparation
1. **Create Seurat Objects**: Prepare two Seurat objects, one for control biopsies (`Control`) and one for wounded biopsies (`Wounded`).
2. **Assign Metadata**: Add a grouping factor to the Seurat objects before integration.

```R
Control$condition <- "Control"
Wounded$condition <- "Wounded"
```

3. **Perform Integration**: Follow the standard Seurat workflow to integrate the datasets, resulting in a combined Seurat object named `Cell.integrated`. Use `"condition"` as the grouping factor in this example.

### Example Analyses

#### Calculating Conserved Expressed Genes
To generate a list of conserved expressed genes across groups:

```R
conserved.genes <- scCEGs(seu_obj = Cell.integrated, grouping.var = "condition")
View(conserved.genes)
```

#### Calculating Differentially Expressed Genes (p-value < 0.05)
To generate a list of DEGs with a p-value < 0.05:

```R
DEGs.pos <- scDEGs(seu_obj = Cell.integrated, ident.1 = "Wounded", ident.2 = "Control", 
                   group.by = "condition", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
View(DEGs.pos)
```

#### Calculating All Differentially Expressed Genes
To generate a list of all DEGs, including those with a p-value ≥ 0.05:

```R
DEGs.all <- scDEGs(seu_obj = Cell.integrated, ident.1 = "Wounded", ident.2 = "Control", 
                   group.by = "condition", min.pct = 0.25, logfc.threshold = 0.25, only.pos = FALSE)
View(DEGs.all)
```

## Notes
- Ensure that the Seurat object provided to `scCEGs()` and `scDEGs()` has undergone integration and contains the specified grouping variable.
- The package automatically excludes clusters with fewer than three cells during analysis, enhancing the reliability of the results.
- For additional support or to report issues, please visit the [GitHub repository](https://github.com/LingzhangMeng/scTidyGene).