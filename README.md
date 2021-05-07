# zinbwave
Zero-inflated Negative Binomial based Wanted Variation Extraction (ZINB-WaVE)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![BioC release](http://www.bioconductor.org/shields/build/release/bioc/zinbwave.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/zinbwave)
[![BioC devel](http://www.bioconductor.org/shields/build/release/bioc/zinbwave.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/zinbwave)
[![R-CMD-check](https://github.com/drisso/zinbwave/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/drisso/zinbwave/actions)
[![Coverage](https://codecov.io/gh/drisso/zinbwave/branch/master/graph/badge.svg)](https://codecov.io/gh/drisso/zinbwave)

This package implements a zero-inflated negative binomial model for single-cell RNA-seq data, with latent factors.

The model is described in details in the paper:

[D. Risso, F. Perraudeau, S. Gribkova, S. Dudoit and JP. Vert (2018).
A general and flexible method for signal extraction from single-cell RNA-seq data. Nature Communications.](https://www.nature.com/articles/s41467-017-02554-5)

## Installation

Since Bioconductor 3.7 the new recommended way to install Bioconductor packages is via the BiocManager package, available on CRAN:

```{r}
install.packages("BiocManager")
BiocManager::install("zinbwave")
```

Note that `zinbwave` requires `R (>=3.4)` and `Bioconductor (>=3.6)`.

In virtually all cases, installing from Bioconductor is recommended. However, if you want to install the development version of `zinbwave` from GitHub, you can do so with the following.

```{r}
library(devtools)
install_github("drisso/zinbwave")
```
