# Description
The R package scTidyGene is developed for geneating conserved expressed genes and differentially expressed genes in scRNA-seq analysis in 
a loop manner, which means automatically calculate and generate list(s) containing conserved expressed genes and/or differentially expressed genes across all clusters, 
. During the calculation and generation the list, this package could deteck and avoid those clusters with less than 3 cells.

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
