# Description
The R package is developed for geneating conserved expressed genes and differentially expressed genes in scRNA-seq analysis in 
a loop manner, which means automatically calculate and generate a list containing conserved expressed genes across all clusters, </br>
and could automatically calculate and generate differentially expressed genes accross all clusters. During the calculation </br>
and generation the list, this package could deteck and avoid those clusters with less than 3 cells.

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
