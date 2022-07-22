# CProtMEDIAS
---
DNA/RNA/Protein sequences general analysis workflow after multiple sequences alignment

A sequence analysis system including sequence digitization, dimensionality reduction, 
specific site search, pseudotime analysis, network construction and result visualization.

### Install
---
if (!require(devtools)){install.packages("devtools")}

devtools::install_github("Busydog1990/CProtMEDIAS")

### vignette
---
vignette("genepro")

If you have any questions, contact me <zzhe@webmail.hzau.edu.cn>

### Note

As from R 4.2.0, conditions of length greater than one are an error in “if” (warnings in previous R version). 
Run monocle_workflow function with R version above 4.2.0 will report an error.
I will fix it in the upcoming version.
