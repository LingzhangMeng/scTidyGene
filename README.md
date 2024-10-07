# Description
The R package scTidyGene is developed for generating list(s) of Conserved Expressed Genes and Differentially Expressed Genes during scRNA-seq analysis in 
a loop manner, meaning automatically calculates and generates list(s) containing Conserved Expressed Genes and/or Differentially Expressed Genes across all cell clusters (cell types). During the calculation and generation of these gene list(s), this package could deteck and avoid those clusters with less than 3 cells, reporting automatically.

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
Function 1: scCEGs()  </br></br>
 to calculate and generate a list for conserved expressed genes</br></br>
scCEGs(seu_obj, grouping.var) </br></br>
seu_obj: seurat object including at least two groups by performing integration.</br>
grouping.var: grouping factor </br></br>

Function : scDEGs()  </br></br>
 to calculate and generate a list for differentially expressed genes</br> </br>
scDEGs(seu_obj, ident.1, ident.2, group.by, min.pct, logfc.threshold = 0.25, only.pos = FALSE) </br>
seu_obj: Character. Seurat object including at least two groups by performing integration. </br>
ident.1: Character. Name of the first group, numerator in the calculation </br>
ident.2: Character. Name of the second group, denominator in the calculation </br>
group.by: Character. Grouping factor </br>
min.pct: Number. Genes expressed in minimal percentage of cells. 0.25 was set as default in the package,  </br>
but it could be customized depends user(s).   </br>
logfc.threshold: number: number, fold change of gene expression level in minimal percentage of cells.  </br>
0.25 was set as default in the package, but it could be customized depends user(s). </br>
only.pos: TRUE/FALSE. Set it to TRUE to select differentially expressed gene list with p value < 0.05; </br>
Set it to FALSE to select all gene list even those with p value >= 0.05. "only.pos = TRUE" was set default in this package.</br>
# Examples
Data Preparaton </br>
Firstly, I have created a seurat object with control biopsies (named Control) and a seurat object with Wounded biopsies (named Wounded).</br>
Secondly, I added assigned new information before integration with the below scripts: </br>
Control$condition <- "Control"     </br>
Wounded$condition <- "Wounded"     </br>
* Please note that I used "condition" as grouping factor in this tutorial;  </br>
* After integration, I generated a seurat object named "Cell.integrated". </br> 
Thirdly, perform integration following standard workflow.</br>
Then </br> </br>
calculate conserved expressed genes and generate a list </br></br>

conserved.gene <- scCEGs(Cell.integrated, grouping.var = "condition") </br></br>
check the gene list with script below: </br>
View(conserved.gene) </br></br>

calculate differentially expressed genes and generate a list including those with p value < 0.05 </br></br>

DEGs.pos <- scDEGs(Cell.integrated, ident.1 = "Wounded", ident.2 = "Control", group.by = "condition", </br>
min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE) </br></br>
check the gene list with script below: </br>
View(DEGs.pos) </br></br>

calculate differentially expressed genes and generate a list including all gene even those with p value >= 0.05 </br></br>

DEGs.all <- scDEGs(Cell.integrated, ident.1 = "Wounded", ident.2 = "Control", group.by = "condition", </br>
min.pct = 0.25, logfc.threshold = 0.25, only.pos = FALSE) </br></br>
check the gene list with script below: </br>
View(DEGs.all) </br></br>
















