# ShrinkNet

This R package implements the method described in

Leday, G. G. R., de Gunst, M. C. M., Kpogbezan, G. B., van der Vaart, A. W., van Wieringen, W. N., and van de Wiel, M. A. (2015).
[Gene network reconstruction using global-local shrinkage priors](http://arxiv.org/abs/1510.03771). *Submitted*.

## Description

**ShrinkNet** enables the reconstruction of an undirected network from high-throughput molecular data. Although it was primarily developed to analyse mRNA expression data, the method is general and can be applied to any data set for which it is reasonable to assume a multivariate Gaussian model. In genomics, this typically includes molecular data generated from a microarray technologies. Hence, the software can in principle be used to analyse protein expression data (e.g. as produced by reverse-phase protein arrays), microRNA and metabolomic data (among others).

**ShrinkNet** aims to be computationally efficient. Core functions are implemented in C++ using the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) software packages, and SVD decompositions are employed to speed up the variational algorithm. Furthermore, the package has been designed so the most computationally intensive steps can be parallelized.

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

## Example 1

- **Data:** gene expression data (B-lymphocyte cells) analysed by Mohammadi and Wit (2015) 
- **Samples:** 60 unrelated individuals
- **Variables:** 100 most variables probes/genes

```R
# Load library
library(BDgraph)

# Load the gene expression data
data(geneExpression)

# Center and scale the data
mytX <- t(scale(geneExpression))

# Run ShrinkNet
res <- ShrinkNet(tX=mytX)
```


## References

Mohammadi, A. and Wit, E. C. (2015). Bayesian structure learning in sparse Gaussian graphical models. *Bayesian Anal*. **10** 109-138.

