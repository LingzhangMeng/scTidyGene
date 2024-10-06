# Description
The R package scTidyGene is developed for generating Conserved Expressed Genes and Differentially Expressed Genes in scRNA-seq analysis in 
a loop manner, which means automatically calculate and generate list(s) containing Conserved Expressed Genes and/or Differentially Expressed Genes across all cell clusters (cell types). During the calculation and generation of these gene list(s), this package could deteck and avoid those clusters with less than 3 cells, reporting automatically.

# Installation
library(devtools) </br>
devtools::install_github("LingzhangMeng/scTidyGene")

# Dependencies
library(Seurat) </br>
library(progress) </br>
library(dplyr) </br>
library(tidyverse) </br>
library(data.table) </br>

# User Guide
