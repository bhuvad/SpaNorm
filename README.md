# SpaNorm - Spatially-aware normalisation for spatial transcriptomics data
<!-- badges: start -->
[![R-CMD-check](https://github.com/bhuvad/spaNorm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bhuvad/spaNorm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

SpaNorm is a spatially aware library size normalisation method that removes library size effects, while retaining biology. Library sizes need to be removed from molecular datasets to allow comparisons across observations, in this case, across space. Bhuva et al. (2024) and Atta et al. (2023) have shown that standard single-cell inspired library size normalisation approaches are not appropriate for spatial molecular datasets as they often remove biological signals while doing so. This is because library size confounds biology in spatial molecular data.

![_The SpaNorm workflow: SpaNorm takes the gene expression data and spatial coordinates as inputs. Using a gene-wise model (e.g., Negative Binomial (NB)), SpaNorm decomposes spatially-smooth variation into those unrelated to library size (LS), representing the underlying true biology and those related to library size. The adjusted data is then produced by keeping only the variation unrelated to library size._](vignettes/SpaNormWorkflow.png)

SpaNorm uses a unique approach to spatially constraint modelling approach to model gene expression (e.g., counts) and remove library size effects, while retaining biology. It achieves this through three key innovations:

1. Computing spatially smooth functions (using thin plate splines) to represent the gene- and location-/cell-/spot- specific size factors.
1. Optmial decomposition of spatial variation into spatially smooth library size associated (technical) and library size independent (biology) variation using generalized linear models (GLMs).
1. Adjustment of data using percentile adjusted counts (PAC) (Salim et al., 2022), as well as other adjustment approaches (e.g., Pearson).

## Installation

SpaNorm can be installed from Bioconductor directly as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SpaNorm")
```
