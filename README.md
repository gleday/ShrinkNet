# ShrinkNet

This R package implements the method described in

Leday, G. G. R., de Gunst, M. C. M., Kpogbezan, G. B., van der Vaart, A. W., van Wieringen, W. N., and van de Wiel, M. A. (2015).
[Gene network reconstruction using global-local shrinkage priors](http://arxiv.org/abs/1510.03771). *Submitted*.

## Description

**ShrinkNet** enables the reconstruction of an undirected network from high-throughput molecular data. Although it was primarily developed to analyse mRNA expression data, the method is general and can be applied to any data set for which it is reasonable to assume a multivariate Gaussian model. In genomics, this typically includes molecular data generated from a microarray technologies. Hence, the software can in principle be used to analyse protein expression data (e.g. as produced by reverse-phase protein arrays), microRNA and metabolomic data (among others).

**ShrinkNet** aims to be computationally efficient. Core functions are implemented in C++ using the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) software packages, and SVD decompositions are employed to speed up the variational algorithm. Furthermore, the package has been designed so the most computationally intensive steps can be parallelized.

## Installation

If you wish to install the package from R:

1) Install and load the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):
```R
install.packages("devtools")
library(devtools)
```
2) Install and load **ShrinkNet**:
```R
install_github("gleday/ShrinkNet")
library(ShrinkNet)
```

## Example 1

- **Data:** gene expression data (B-lymphocyte cells) analysed by Mohammadi and Wit (2015) 
- **Samples:** 60 unrelated individuals
- **Variables:** 100 most variables probes/genes

R commands:

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

Output:

```
STEP 0: SVD computations... DONE
STEP 1: Variational algorithm...
iteration 1
iteration 2
iteration 3
iteration 4
iteration 5
...
iteration 74
DONE
STEP 2: Calculate summary statistics from posteriors... DONE
STEP 3: Estimate p0... DONE
STEP 4: Edge selection... DONE

prior null probability p0 = 0.70242 
73 selected edges out of 4950 (1.47%) using blfdr = 0.1

Time (H:MM:SS): 0:00:16
```

## Example 2

- **Data:** Protein expression data from [TCPA](http://app1.bioinformatics.mdanderson.org/tcpa/_design/basic/index.html) - Ovarian serous cystadenocarcinoma (OV)
- **Samples:** 412 tumor samples
- **Variables:** 190 proteins

R commands:

```R
# Read data
tcpaOV <- read.table(file="TCGA-OV-L3-S35.csv", sep=",", header=TRUE)

# Remove annotations
datamatOV <- data.matrix(tcpaOV[,-c(1,2)])
rownames(datamatOV) <- tcpaOV$TCGA_patient_barcode

# Center and scale the data
mytX <- t(scale(datamatOV, center = TRUE, scale = TRUE))

# Run ShrinkNet
res <- ShrinkNet(tX=mytX, methodp0="sampling", nsamp=1000, ncpus=8)
```

Output:

```
STEP 0: SVD computations... 
R Version:  R version 3.1.1 (2014-07-10) 

snowfall 1.84-6.1 initialized (using snow 0.3-13): parallel execution on 8 CPUs.

Library ShrinkNet loaded.
Library ShrinkNet loaded in cluster.


Stopping cluster

DONE
STEP 1: Variational algorithm...
iteration 1
iteration 2
iteration 3
iteration 4
iteration 5
iteration 6
iteration 7
DONE
STEP 2: Calculate summary statistics from posteriors... 
snowfall 1.84-6.1 initialized (using snow 0.3-13): parallel execution on 8 CPUs.

Library ShrinkNet loaded.
Library ShrinkNet loaded in cluster.


Stopping cluster

DONE
STEP 3: Estimate p0... 
snowfall 1.84-6.1 initialized (using snow 0.3-13): parallel execution on 8 CPUs.

Library ShrinkNet loaded.
Library ShrinkNet loaded in cluster.


Stopping cluster

DONE
STEP 4: Edge selection... DONE

prior null probability p0 = 0.896 
923 selected edges out of 17955 (5.14%) using blfdr = 0.1

Time (H:MM:SS): 0:00:54
```

## References

Mohammadi, A. and Wit, E. C. (2015). Bayesian structure learning in sparse Gaussian graphical models. *Bayesian Anal*. **10** 109-138.

