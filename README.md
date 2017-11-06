[![Travis-CI Build Status](https://travis-ci.org/arendsee/phylostratr.svg?branch=master)](https://travis-ci.org/arendsee/phylostratr)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/phylostratr/master.svg)](https://codecov.io/github/arendsee/phylostratr?branch=master)

# phylostratr

## Building a database

The goal of phylostratigraphy is to pinpoint the origin of a gene. We can do
this by searching for homologs within increasingly broad clades. The highest
clade that contains all homologs of the gene is the gene's phylostratum. For
example, if a gene has homologs across members of Brassicaceae, but none
outside, then we say it is a Brassicaceae-specific gene and a member of the
Brassicaceae phylostratum.

I will first dicuss ideal phylostratigraphy, where homology can be inferred, if
present, with high certainty across all species of an accurate species tree. 

Issues in the ideal case:

 1. chimerism
 2. deletion
 3. lateral transfer
 4. convergent evolution
 5. sparse trees

But in reality, there are many additional complexities:

 1. imperfect homology inference
 2. incomplete sequence data
 3. random uncertainty in gene annotation
 4. systematic differences in gene annotation 
 5. uncertainty in tree topology

The NCBI common tree is very low resolution, with many multifurcating nodes.

Defining a search space.

While it may be tempting to just BLAST your sequences against all of nr, this
approach is not suitable for phylostratigraphy.
