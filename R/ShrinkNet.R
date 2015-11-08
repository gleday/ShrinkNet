#' Gene network reconstruction using global-local shrinkage priors
#'
#' @param tX p by n data matrix of gene expression measurements (centered and scaled)
#' @param globalShrink integer. Either 1 or 2. See Details.
#' @param methodp0 character. Either "exact" or "sampling". See Details.
#' @param nsamp integer. Number of randomly selected edges to estimate p0. Only when methodp0="sampling". See Details.
#' @param ncpus integer. Number of cpu cores to be used.
#' @param blfdr numeric. Bayesian analogue of the local false discovery rate used for edge selection. Value should be between 0 and 1. Default is 0.1.
#' @param maxiter integer. Maximum number of iterations for the variational algorithm. Default is 100.
#' @param tol numeric. Represents the maximum relative convergence tolerance over the p variational lower bounds. Default is 0.001.
#' @param verbose logical. Should information on progress be printed?
#' @param standardize logical. Should the data be standardized? Default is TRUE.
#' @details
#' The function enables the reconstruction of an undirected network given data. 
#' Although the tool was primarily developed to analyse mRNA expression data,
#' the method is general and can be applied to any data set for which it is reasonable to
#' assume a multivariate Gaussian model. This typically include molecular data generated from
#' a microarray-based technology such as protein expression data (e.g. as produced by reverse-phase protein arrays)
#' and microRNA data.
#' 
#' \code{ShrinkNet} aims to be computationally efficient. Core functions are implemented in C++ using
#' the Rcpp and RcppArmadillo software packages and SVD decompositions are employed to speed up the
#' variational algorithm. Furthermore, functions in the package has been designed so that the most
#' computationally intensive steps can be parallelized.
#' 
#' If \code{globalShrink}=1 then empirical Bayes for the global shrinkage prior is carried out using
#' fixed-point iterations as in Valpola and Honkela (2006). Otherwise, if \code{globalShrink}=2,
#' the approximate analytical solution of Leday et al (2015) is used.
#'
#' If \code{methodp0}="exact", all P=0.5*p*(p-1) edges are investigated whereas if \code{methodp0}="sampling"
#' a random subset of size \code{nsamp} is selected to estimate p0. When p<=100 we suggest using \code{methodp0}="exact" and
#' \code{methodp0}="sampling" otherwise with \code{nsamp} being greater than 1000.
#' 
#' @return An object of class \code{\link{ShrinkNet-class}}
#' 
#' @author Gwenael G.R. Leday <gwenael.leday (at) mrc-bsu.cam.ac.uk>
#' 
#' @references Leday, G.G.R., de Gunst, M.C.M., Kpogbezan, G.B., van der Vaart, A.W., van Wieringen, W.N., and
#' van de Wiel, M.A. (2015). Gene network reconstruction using global-local shrinkage priors. Submitted.
#' 
#' @export
ShrinkNet <- function(tX, globalShrink=1, methodp0="exact", nsamp=1000, ncpus=1, blfdr=0.1, maxiter=100, tol=0.001, verbose=TRUE, standardize=TRUE){

  ##### Input checks
  if(is.matrix(tX)){
    if(any(is.na(tX))){
      stop("Missing values are not allowed")
    }
  }else{
    stop("tX is not a matrix")
  }
  if(is.numeric(globalShrink)){
    globalShrink <- as.integer(globalShrink)
    if(!globalShrink%in%c(1,2)){
      stop("globalShrink should be equal to 1 or 2")
    }
  }else{
    stop("globalShrink is not a numeric")
  }
  if(is.character(methodp0)){
    if(!methodp0%in%c("exact","sampling")){
      stop("methodp0 should be equal to 'exact' or 'sampling'.")
    }
    if(methodp0=="exact" & min(dim(tX))>100){
      methodp0 <- "sampling"
      warning("min(n,p)>100 so methodp0 has been changed to 'sampling'")
    }
  }else{
    stop("methodp0 is not a character")
  }
  edgeTot <- 0.5*nrow(tX)*(nrow(tX)-1)
  if(is.numeric(nsamp)){
    nsamp <- as.integer(nsamp)
    if(nsamp>edgeTot){
      nsamp <- edgeTot
    }
    if(nsamp<1000){
      warning("nsamp (<1000) may be too low to obtain a reasonable estimate of p0")
    }
  }else{
    stop("nsamp is not a numeric")
  }
  if(is.numeric(maxiter)){
    maxiter <- as.integer(maxiter)
  }else{
    stop("maxiter is not a numeric")
  }
  if(is.numeric(blfdr)){
    if((blfdr<=0)|(blfdr>=1)){
      stop("blfdr should be between 0 and 1")
    }
  }else{
    stop("blfdr is not a numeric")
  }
  if(is.numeric(ncpus)){
    ncpus <- as.integer(ncpus)
    if(ncpus>1){
      if(max(dim(tX))<100){
        ncpus <- 1
        warning("max(n,p)<100: no parallel computations")
      }
    }
  }else{
    stop("ncpus is not a numeric")
  }
  if(!is.logical(verbose)){
    stop("verbose is not a logical")
  }
  if(!is.logical(standardize)){
    stop("standardize is not a logical")
  }else{
    if(standardize){
      tX <- t(scale(t(tX), center = TRUE, scale = TRUE))
    }else{
      warning("Input data have not been standardized")
    }
  }
  
  tps <- proc.time()

  ##### Initialization
  aRand <- 0.001
  bRand <- 0.001

  ##### Data preparation
  if(verbose){
    cat("\n")
    cat("STEP 0: SVD computations... ")
  }
  if(ncpus==1){
    allSVDs <- sapply(1:nrow(tX), getSVD, tX=tX, simplify=FALSE)
  }else{
    cat("\n")
    snowfall::sfInit(parallel=TRUE, cpus=ncpus)
    snowfall::sfLibrary(ShrinkNet)
    allSVDs <- snowfall::sfSapply(1:nrow(tX), getSVD, tX=tX, simplify=FALSE)
    snowfall::sfRemoveAll()
    snowfall::sfStop()
  }
  if(verbose) cat("DONE\n")
  
  ##### Algo
  if(verbose) cat("STEP 1: Variational algorithm...\n")
  eb <- HiddenVarAlgo(SVDs=allSVDs, aRand=aRand, bRand=bRand, maxiter=maxiter, globalShrink=globalShrink, tol=tol, verbose=verbose)

  ##### Calculate summary statistics from posteriors
  if(verbose) cat("STEP 2: Calculate summary statistics from posteriors... ")
  if(ncpus==1){
    matThres <- sapply(1:nrow(tX), HiddenVarRidgeiGetKappa, SVDs=allSVDs, aRand=eb$parTau[nrow(eb$parTau),1], bRand=eb$parTau[nrow(eb$parTau),2], bRandStarInit=eb$allbRandStar, dSigmaStarInit=eb$alldSigmaStar, simplify=TRUE)
  }else{
    cat("\n")
    snowfall::sfInit(parallel=TRUE, cpus=ncpus)
    snowfall::sfLibrary(ShrinkNet)
    matThres <- snowfall::sfSapply(1:nrow(tX), HiddenVarRidgeiGetKappa, SVDs=allSVDs, aRand=eb$parTau[nrow(eb$parTau),1], bRand=eb$parTau[nrow(eb$parTau),2], bRandStarInit=eb$allbRandStar, dSigmaStarInit=eb$alldSigmaStar, simplify=TRUE)
    snowfall::sfRemoveAll()
    snowfall::sfStop()
  }
  matThres <- (matThres + t(matThres))/2
  if(verbose) cat("DONE\n")
  
  ##### Estimate p0
  if(verbose) cat("STEP 3: Estimate p0... ")
  if(methodp0=="exact"){
    p0 <- HiddenEstimatep0(themat=matThres, tX=tX)
    #p0 <- resp0$p0
  }else{
    mat <- matThres
    mat[upper.tri(mat)] <- 0
    idx <- which(mat!=0, arr.ind=TRUE)
    idx <- idx[sample(nrow(idx), nsamp),]
    if(ncpus==1){
      allLogBFs <- t(apply(idx, 1, .edgeBFprime, themat=matThres, tX=tX))
    }else{
      cat("\n")
      snowfall::sfInit(parallel=TRUE, cpus=ncpus)
      snowfall::sfLibrary(ShrinkNet)
      allLogBFs <- t(snowfall::sfApply(idx, 1, .edgeBFprime, themat=matThres, tX=tX))
      snowfall::sfRemoveAll()
      snowfall::sfStop()
    }
    p0 <- 1-mean(allLogBFs>0)
  }
  if(verbose) cat("DONE\n")

  ##### Edge selection using Bayesian local false discovery rate
  if(verbose) cat("STEP 4: Edge selection... ")
  selGraph <- HiddenEdgeSelection(themat=matThres, tX=tX, p0=p0, lfdrcut=blfdr)
  nbedge <- sum(selGraph)/2
  if(verbose){
    cat("DONE\n\n")
    cat("prior null probability p0 =", round(p0,5), "\n")
    cat("", nbedge, " selected edges out of ", edgeTot, " (",round(100*nbedge/edgeTot, 2),"%)", " using blfdr = ", blfdr, sep="")
  }
  tps2 <- proc.time() - tps
  if(verbose){
    cat("\n\n")
    cat("Time (H:MM:SS):", .convertToTime(tps2[3]))
    cat("\n\n")
  }
  
  ## Output
  myigraph <- igraph::graph.adjacency(selGraph, mode = "undirected")
  if(!is.null(rownames(tX))){
    myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=rownames(tX))
  }
  out <- new("ShrinkNet",
             graph = myigraph,
             kappa = matThres,
             p0 = p0,
             globalPrior = eb$parTau,
             allmargs = eb$allmargs,
             time = tps2[3])
  
  return(out)
}
