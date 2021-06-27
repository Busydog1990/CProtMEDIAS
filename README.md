# genepro

DNA/RNA/Protein sequences general analysis workflow after multiple sequences alignment

A sequence analysis system including sequence digitization, dimensionality reduction, 
specific site search, pseudotime analysis, network construction and result visualization.

### Install

if (!require(BiocManager)){install.packages("BiocManager")}

required_packages <- c("Biostrings", "Seurat", "cowplot", "ggplot2", "ggrepel","ggsci", "ggseqlogo",
                       "igraph", "monocle", "reshape2", "VGAM", "BiocGenerics", "xtable", "SeuratObject","devtools")

BiocManager::install(required_packages[!required_packages %in% rownames(installed.packages())])

devtools::install_github("Busydog1990/genepro")

### vignette

vignette("genepro")

If you have any questions, contact me <zzhe@webmail.hzau.edu.cn>
