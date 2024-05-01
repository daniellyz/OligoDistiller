
# Installation from Github in Rstudio

```R
# Install BiocManager if it has not been installed previously:
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# Install OligoDistiller:
BiocManager::install("daniellyz/OligoDistiller")
```
# Webtool

https://valkenborg-lab.shinyapps.io/OligoDistiller/

# Examples for both webtool and R package:

## Tutorial 1: Oligonucleotide Impurity profiling from a crude sample

The tutorial is available at: https://github.com/daniellyz/mergeion.github.io/blob/7e9fa0eb2c34060dd6eb622053d2c53d04f2b021/Tutorial%20MS1.pdf

## Tutorial 2: ESI-MS/MS for sequence confirmation of a synthetic oligonucleotide

The tutorial is available at:
https://github.com/daniellyz/mergeion.github.io/blob/7e9fa0eb2c34060dd6eb622053d2c53d04f2b021/Tutorial%20ESI_MS2.pdf

## Tutorial 3: MALDI-FT-MS/MS for sequence confirmation of a simple oligonucleotide

The tutorial is available at:
https://github.com/daniellyz/mergeion.github.io/blob/7e9fa0eb2c34060dd6eb622053d2c53d04f2b021/Tutorial%20MALDI-MS2.pdf
