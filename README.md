# ShrinkNet

This R package implements the methods described in

Leday, G. G. R., de Gunst, M. C. M., Kpogbezan, G. B., van der Vaart, A. W., van Wieringen, W. N., and van de Wiel, M. A. (2016).
[Gene network reconstruction using global-local shrinkage priors](http://www.e-publications.org/ims/submission/AOAS/user/submissionFile/25409?confirm=c2e87384). *The Annals of Applied Statistics (to appear)*.

## Description

ShrinkNet enables the reconstruction of an undirected network from high-throughput molecular data that characterizes the *conditional independence* structure between the molecular variables.

Although ShrinkNet was primarily developed to analyse mRNA expression data, the method is general and can be applied to any data for which it is reasonable to assume a multivariate Gaussian model. In genomics, this may include protein ([Abkani et al., 2014](http://dx.doi.org/10.1038/ncomms4887)), microRNA ([Stingo et al., 2010](http://projecteuclid.org/euclid.aoas/1294167808)) and metabolomic data ([Krumsiek et al., 2011](http://www.biomedcentral.com/1752-0509/5/21)) for example.

ShrinkNet aims to be computationally efficient. Core functions are implemented in C++ using the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html) software packages, and SVD decompositions are employed to speed up the variational algorithm.


## Installation

If you wish to install **ShrinkNet** from R (using package [devtools](https://cran.r-project.org/web/packages/devtools/index.html)):

```R
install.packages("devtools")
library(devtools)
install_github("gleday/ShrinkNet")
library(ShrinkNet)
```

**Remark:** For Mac OS X, it may be necessary to install a recent GNU Fortran compiler beforehand. See [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/) for more explanations.

## Usage

Given a data matrix *mydata*, where rows represent molecular variables (e.g. genes) and columns represent samples, ShrinkNet is run by the following R command:
```R
myobject <- ShrinkNet(tX=mydata)
```

The R function ShrinkNet() returns an (S4) object of class "ShrinkNet" that is associated with the following functions:

```R
# Print object information
print(myobject)

# Summary on inferred network
summary(myobject)

# Extract adjacency matrix corresponding to the inferred network
adjacency(myobject)

# Extract matrix containing scores on edges (and used to rank them)
score(myobject)

# Extract matrix containing the posterior means of regression coefficients in the SEM
coef(myobject)

# Plot convergence of the variational algorithm
plotML(myobject)

# Plot density of the global shrinkage prior
plotPrior(myobject)

# Plot the inferred graph
plotGraph(myobject)

# Extract a table of the selected edges
listEdges(myobject)

# Extract a table of the 'nb' top-ranked edges
topEdges(myobject, nb=10)

# Extract a table of the 'nb' top-ranked nodes
topDegree(myobject, nb=10)
```


## Example 1

- **Data:** gene expression data (B-lymphocyte cells) analysed by [Mohammadi and Wit (2015)](http://projecteuclid.org/euclid.ba/1422468425)
- **Samples:** 60 unrelated individuals
- **Variables:** 100 most variables probes/genes

R commands:

```R
# Load library
library(BDgraph)

# Load the gene expression data
data(geneExpression)

# Run ShrinkNet
res <- ShrinkNet(tX=t(geneExpression))
```

R console:

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
78 selected edges out of 4950 (1.58%) using blfdr = 0.1

Time (H:MM:SS): 0:00:17
```

## Example 2

- **Data:** gene expression data from R package [care](https://cran.r-project.org/web/packages/care/index.html)
- **Samples:** 30 human brain samples
- **Variables:** 403 genes

R commands:

```R
# Load library
library(care)

# Load the gene expression data
data(lu2004)

# Run ShrinkNet
res <- ShrinkNet(tX=t(lu2004$x), nsamp0=10000)
```

R console:

```
STEP 0: SVD computations... DONE
STEP 1: Variational algorithm...
iteration 1
iteration 2
iteration 3
iteration 4
iteration 5
iteration 6
iteration 7
iteration 8
iteration 9
iteration 10
DONE
STEP 2: Calculate summary statistics from posteriors... DONE
STEP 3: Estimate p0... DONE
STEP 4: Edge selection... DONE

prior null probability p0 = 0.75775 
916 selected edges out of 81003 (1.13%) using blfdr = 0.1

Time (H:MM:SS): 0:00:58
```


## Example 3

- **Data:** gene expression data from R package [GeneNet](https://cran.r-project.org/web/packages/GeneNet/index.html)
- **Time points:** 9
- **Variables:** 102 genes

R commands:

```R
# Load library
library(GeneNet)

# Load the gene expression data
data(ecoli)

# Run ShrinkNet
res <- ShrinkNet(tX=t(ecoli))
```

R console:

```
STEP 0: SVD computations... DONE
STEP 1: Variational algorithm...
iteration 1
iteration 2
iteration 3
iteration 4
iteration 5
iteration 6
iteration 7
iteration 8
iteration 9
DONE
STEP 2: Calculate summary statistics from posteriors... DONE
STEP 3: Estimate p0... DONE
STEP 4: Edge selection... DONE

prior null probability p0 = 0.72996 
94 selected edges out of 5151 (1.82%) using blfdr = 0.1

Time (H:MM:SS): 0:00:02
```


## Example 4

- **Data:** Protein expression data from [TCPA](http://app1.bioinformatics.mdanderson.org/tcpa/_design/basic/index.html) - Ovarian serous cystadenocarcinoma (OV)
- **Samples:** 412 tumor samples
- **Variables:** 190 proteins (as measured by reverse-phase protein arrays)

R commands:

```R
# Read data (we assume the file 'TCGA-OV-L3-S35.csv' has been previously downloaded)
tcpaOV <- read.table(file="TCGA-OV-L3-S35.csv", sep=",", header=TRUE)

# Remove annotations
datamatOV <- data.matrix(tcpaOV[,-c(1,2)])
rownames(datamatOV) <- tcpaOV$TCGA_patient_barcode

# Run ShrinkNet
res <- ShrinkNet(tX=t(datamatOV), nsamp0=10000)
```

R console:

```
STEP 0: SVD computations... DONE
STEP 1: Variational algorithm...
iteration 1
iteration 2
iteration 3
iteration 4
iteration 5
iteration 6
iteration 7
DONE
STEP 2: Calculate summary statistics from posteriors... DONE
STEP 3: Estimate p0... DONE
STEP 4: Edge selection... DONE

prior null probability p0 = 0.8937 
937 selected edges out of 17955 (5.22%) using blfdr = 0.1

Time (H:MM:SS): 0:12:48
```

