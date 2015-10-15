#' Undirected network inference using a Bayesian SEM and a posteriori edge selection
#'
#' @param tX p by n data matrix
#' @param globalShrink either 1 or 2. See Details.
#' @param maxiter integer. Maximum number of iterations for the variational algorithm. Default is 100.
#' @param blfdr Bayesian analogue of the local false discovery rate used for edge selection. Value should be between 0 and 1. Default is 0.1.
#' @param maxedges integer. Maximum number of edges to consider for forward selection. Default is 0.25*p*(p-1). See Details.
#' @details
#' If \code{globalShrink}=1 then empirical Bayes for the global shrinkage prior is carried out using
#' fixed-point iterations as in Valpola and Honkela (2006). Otherwise, if \code{globalShrink}=2,
#' the approximate analytical solution of Leday et al (2015) is used.
#'
#' The edge selection procedure respects the ranking on edges in terms of the order in which
#' they are considered in the forward selection. However, because the Bayes factor generally
#' decreases with the rank, it may be pratical to avoid investigating all P=0.5*p*(p-1)
#' possible edges since a large part of them are not important. The argument \code{maxedges}
#' specifies the maximum number of edges to consider both for calculating p0 and the forward
#' selection procedure. The default is 0.5*P.
#' @return A named list with the following elements:
#'  \item{adjacency}{A sparse matrix containing the adjacency matrix corresponding to the selected graph.}
#'  \item{p0}{A number giving the estimate of the proportion of null hypothesis.}
#'  \item{kappabar}{A matrix containing scores used to rank edges.}
#'  \item{globalPrior}{A matrix containing the (shape and rate) parameters of the global shrinkage prior
#'  (gamma distribution) at each iteration of the variational algorithm.}
#'  \item{allmargs}{A vector containing the values of the variational lower bounds for each regression
#'  equation in the Bayesian SEM at each iteration of the variational algorithm.}
#'  \item{time}{Running time of ShrinkNet.}
#' @author Gwenael G.R. Leday <gwenael.leday (at) mrc-bsu.cam.ac.uk>
#' @export
ShrinkNet <- function(tX, globalShrink=1, maxiter=100, blfdr=0.1, maxedges=NULL){

  ##### Input checks
  if(!is.matrix(tX)){
    stop("tX is not a matrix")
  }
  if(is.numeric(globalShrink)){
    if(!any(globalShrink==c(1,2))){
      stop("globalShrink should be equal to 1 or 2")
    }
  }else{
    stop("globalShrink is not a numeric")
  }
  if(!is.numeric(maxiter)){
    stop("maxiter is not a numeric")
  }
  maxiter <- round(maxiter)
  if(!is.null(maxedges)){
    if(!is.numeric(maxedges)){
      stop("maxedges is not a numeric")
    }
    maxedges <- round(maxedges)
  }else{
    maxedges <- round(0.5*nrow(tX)*(nrow(tX)-1)/2)
  }
  if(is.numeric(blfdr)){
    if((blfdr<=0)|(blfdr>=1)){
      stop("blfdr should be between 0 and 1")
    }
  }else{
    stop("blfdr is not a numeric")
  }

  tps <- proc.time()

  ##### Initialization
  aRand <- 0.001
  bRand <- 0.001

  ##### Data preparation
  cat("\n")
  cat("STEP 0: SVD computations... ")
  allSVDs <- sapply(1:nrow(tX), .getSVD, tX=tX, simplify=FALSE)
  cat("DONE\n")

  ##### Algo
  cat("STEP 1: Variational algorithm...\n")
  eb <- varAlgo(SVDs=allSVDs, aRand=aRand, bRand=bRand, maxiter=maxiter, globalShrink=globalShrink)

  ##### Estimate p0
  cat("STEP 2: Estimate p0... ")
  p0 <- estimatep0(themat=eb$matThres, tX=tX, maxedges=maxedges)
  cat("DONE\n")

  ##### Edge selection using Bayesian local false discovery rate
  cat("STEP 3: Edge selection... ")
  selGraph <- edgeSelection(themat=eb$matThres, tX=tX, p0=p0, maxedges=maxedges, lfdrcut=blfdr)
  cat("DONE\n\n")

  cat("prior null probability p0 = ", round(p0,5), "\n")
  nbedge <- sum(selGraph)/2
  edgeTot <- 0.5*nrow(eb$matThres)*(nrow(eb$matThres)-1)
  cat("", nbedge, " selected edges out of ", edgeTot, " (",round(100*nbedge/edgeTot, 2),"%)", "", sep="")
  tps2 <- proc.time() - tps
  cat("\n\n")
  cat("Time (steps 0 to 3):", .convertToTime(tps2[3]))
  #print(tps2)
  cat("\n\n")

  ## Output
  out <- list(
    "adjacency" = Matrix(selGraph, sparse=TRUE),
    "p0"= p0,
    "kappabar" = eb$matThres,
    "globalPrior" = eb$parTau,
    "allmargs" = eb$allmargs,
    "time" = tps2)
  return(out)
}
