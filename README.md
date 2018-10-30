[![Travis-CI Build Status](https://travis-ci.org/arendsee/phylostratr.svg?branch=master)](https://travis-ci.org/arendsee/phylostratr)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/phylostratr/master.svg)](https://codecov.io/github/arendsee/phylostratr?branch=master)
[![DOI](https://zenodo.org/badge/109036472.svg)](https://zenodo.org/badge/latestdoi/109036472)

# phylostratr

Predict and explore the age of genes using phylostratigraphic methods.

![Phylostratr Workflow](./README-fig1.png)

## Installation

You can install from GitHub with:

```{r github-installation, eval=FALSE}
library(devtools)
install_github("arendsee/phylostratr")
```

The above command currently fails with a cryptic 404 error. I haven't tracked
down the bug yet, but below is an ugly work around:

```R
# from bioconductor
BiocManager::install('Biostrings')

# from github
devtools::install_github('ropensci/taxizedb')

# from CRAN (feel free to ignore any packages that are already installed)
install.packages('ape')
install.packages('curl')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('scales')
install.packages('glue')
install.packages('gridExtra')
install.packages('magrittr')
install.packages('methods')
install.packages('purrr')
install.packages('readr')
install.packages('reshape2')
install.packages('rhmmer')
install.packages('rlang')
install.packages('tibble')
install.packages('tidyr')

# finally install phylostratr itself
install_github("arendsee/phylostratr", dependencies=FALSE)
```

## Dependencies

 * `NCBI BLAST+` - `blastp` (the protein BLAST command) must be in `PATH`. You
   can tell if `blastp` is properly installed by calling `blastp -help` from
   the command line.

## Citation

    Zebulun Arendsee, Jing Li, Urminder Singh, Arun Seetharam, Karin Dorman,
    Eve Syrkin Wurtele (2018) phylostratr: A framework for phylostratigraphy.
    bioRxiv doi: https://doi.org/10.1101/360164

## Funding

This work is funded by the National Science Foundation grant:

[NSF-IOS 1546858](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1546858)
Orphan Genes: An Untapped Genetic Reservoir of Novel Traits
