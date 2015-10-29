# ShrinkNet

This R package implements the method described in

Leday, G. G. R., de Gunst, M. C. M., Kpogbezan, G. B., van der Vaart, A. W., van Wieringen, W. N., and van de Wiel, M. A. (2015).
[Gene network reconstruction using global-local shrinkage priors](http://arxiv.org/abs/1510.03771). *Submitted*.

## Description

**ShrinkNet** enables the reconstruction of an undirected network from high-throughput molecular data. 
Although it was primarily developed to analyse mRNA expression data, the method is general and
can be applied to any data set for which it is reasonable to assume a multivariate Gaussian model.
In genomics, this typically include molecular data generated from a microarray technologies.
Hence, **ShrinkNet** can in principle be used to analyse protein expression data (e.g. as produced by
reverse-phase protein arrays), microRNA and metabolomic data.

**ShrinkNet** aims to be computationally efficient. Core functions are implemented in C++ using
the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) software packages and
SVD decompositions are employed to speed up the variational algorithm.
Furthermore, functions in the package has been designed so that the most
computationally intensive steps can be parallelized.

## Installation

If you wish to install the package from R:

* Install and load the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```R
install.packages("devtools")
library(devtools)
```

* Install and load **ShrinkNet**:

```R
install_github("gleday/ShrinkNet")
library(ShrinkNet)
```


