# zinbwave
Zero-inflated Negative Binomial based Wanted Variation Extraction (ZINB-WaVE)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/drisso/zinbwave.svg?branch=master)](https://travis-ci.org/drisso/zinbwave)
[![Coverage](https://codecov.io/gh/drisso/zinbwave/branch/master/graph/badge.svg)](https://codecov.io/gh/drisso/zinbwave)

This package implements a zero-inflated negative binomial model for single-cell RNA-seq data, with latent factors.

The model is described in details in the paper:

D. Risso, F. Perraudeau, S. Gribkova, S. Dudoit and JP. Vert (2017).
ZINB-WaVE: A general and flexible method for signal extraction from single-cell RNA-seq data. bioRxiv. https://doi.org/10.1101/125112

## Installation

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("zinbwave")
```

Note that `zinbwave` requires `R (>=3.4)` and `Bioconductor (>=3.6)`.

