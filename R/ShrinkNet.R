#' Gene network reconstruction using global-local shrinkage priors
#'
#' @param tX p by n data matrix of gene expression measurements
#' @param globalShrink integer. Either 1 or 2. See Details.
#' @param nsamp0 integer. Number of randomly selected edges to estimate p0. See Details.
#' @param blfdr numeric. Bayesian analogue of the local false discovery rate used for edge selection. Value should be between 0 and 1. Default is 0.1.
#' @param maxNbEdges numeric. Maximum number of edges to select.
#' @param maxiter integer. Maximum number of iterations for the variational algorithm. Default is 100.
#' @param tol numeric. Represents the maximum relative convergence tolerance over the p variational lower bounds. Default is 0.001.
#' @param verbose logical. Should information on progress be printed?
#' @param standardize logical. Should the data be standardized? Default is TRUE.
#' @details
#' 
#' If \code{globalShrink}=1 then empirical Bayes for the global shrinkage prior is carried out using
#' fixed-point iterations as in Valpola and Honkela (2006). Otherwise, if \code{globalShrink}=2,
#' the approximate analytical solution of Leday et al (2015) is used.
#'
#' When \code{nsamp0}=\code{NULL}, the proportion of null hypotheses p0 is estimated using Bayes factors calculated for all P=0.5*p*(p-1) edges (cf Leday et al., 2015).
#' When P is very large, it may me preferable to approximate p0 instead using a random subset of edges.
#' If \code{nsamp0} is an integer, then a random subset of size \code{nsamp0} is selected to estimate p0.
#' The default is \code{nsamp0}=\code{NULL} when p<=100 and \code{nsamp0}=1000 otherwise.
#' 
#' @return An object of class \code{\link{ShrinkNet-class}}
#' 
#' @author Gwenael G.R. Leday <gwenael.leday (at) mrc-bsu.cam.ac.uk>
#' 
#' @references Leday, G.G.R., de Gunst, M.C.M., Kpogbezan, G.B., van der Vaart, A.W., van Wieringen, W.N., and
#' van de Wiel, M.A. (2015). Gene network reconstruction using global-local shrinkage priors. Submitted.
#' 
#' @export
ShrinkNet <- function(tX, globalShrink=1, nsamp0=NULL, blfdr=0.1, maxNbEdges=NULL, maxiter=100, tol=0.001, verbose=TRUE, standardize=TRUE){

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
  edgeTot <- 0.5*nrow(tX)*(nrow(tX)-1)
  if(is.null(nsamp0)){
    if(nrow(tX)>100){
      nsamp0 <- 1000
      warning("p>100 so p0 is estimated by sampling nsamp0=1000 edges")
    }
  }else{
    if(is.numeric(nsamp0)){
      nsamp0 <- as.integer(nsamp0)
      if(nsamp0>edgeTot){
        nsamp0 <- edgeTot
      }
      if(nsamp0<1000){
        warning("nsamp0 (<1000) may be too low to obtain a reasonable estimate of p0")
      }
    }else{
      stop("nsamp0 is not a numeric")
    }
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
  if(is.null(maxNbEdges)){
    maxNbEdges <- 0
  }else{
    if(is.numeric(maxNbEdges)){
      maxNbEdges <- round(maxNbEdges)
      if((maxNbEdges<=0) | (maxNbEdges>edgeTot) ){
        stop(paste("maxNbEdges must take values between 1 and", edgeTot) )
      }
    }else{
      stop("maxNbEdges is not a numeric")
    }
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
  allSVDs <- sapply(1:nrow(tX), getSVD, tX=tX, simplify=FALSE)
  if(verbose) cat("DONE\n")
  
  tps1 <- proc.time() - tps
  
  ##### Algo
  if(verbose) cat("STEP 1: Variational algorithm...\n")
  eb <- HiddenVarAlgo(SVDs=allSVDs, tX=tX, aRand=aRand, bRand=bRand, maxiter=maxiter, globalShrink=globalShrink, tol=tol, verbose=verbose)
  
  tps2 <- proc.time() - tps - tps1

  ##### Calculate summary statistics from posteriors
  if(verbose) cat("STEP 2: Calculate summary statistics from posteriors... ")
  matThres <- sapply(1:nrow(tX), HiddenVarRidgeiGetKappa, SVDs=allSVDs, tX=tX, aRand=eb$parTau[nrow(eb$parTau),1], bRand=eb$parTau[nrow(eb$parTau),2], bRandStarInit=eb$allbRandStar, dSigmaStarInit=eb$alldSigmaStar, simplify=TRUE)
  matThres <- (matThres + t(matThres))/2
  if(verbose) cat("DONE\n")
  
  tps3 <- proc.time() - tps - tps1 - tps2
  
  ##### Estimate p0
  if(verbose) cat("STEP 3: Estimate p0... ")
  if(is.null(nsamp0)){
    p0 <- HiddenEstimatep0(themat=matThres, tX=tX)
  }else{
    mat <- matThres
    mat[upper.tri(mat)] <- 0
    idx <- which(mat!=0, arr.ind=TRUE)
    idx <- idx[sample(nrow(idx), nsamp0),]
    allLogBFs <- t(apply(idx, 1, .edgeBFprime, themat=matThres, tX=tX))
    p0 <- 1-mean(allLogBFs>0)
  }
  if(verbose) cat("DONE\n")
  
  tps4 <- proc.time() - tps - tps1 - tps2 - tps3

  ##### Edge selection using Bayesian local false discovery rate
  if(verbose) cat("STEP 4: Edge selection... ")
  selGraph <- HiddenEdgeSelection(themat=matThres, tX=tX, p0=p0, lfdrcut=blfdr, maxNbEdges=maxNbEdges)
  nbedge <- sum(selGraph)/2
  if(verbose){
    cat("DONE\n\n")
    cat("prior null probability p0 =", round(p0,5), "\n")
    cat("", nbedge, " selected edges out of ", edgeTot, " (",round(100*nbedge/edgeTot, 2),"%)", " using blfdr = ", blfdr, sep="")
  }
  tps5 <- proc.time() - tps - tps1 - tps2 - tps3 - tps4
  tps6 <- proc.time() - tps
  
  ## Time 
  mytime <- data.frame("elapsed"=c(tps1[3], tps2[3], tps3[3], tps4[3], tps5[3], tps6[3]))
  mytime$"H:MM:SS" <- sapply(mytime$elapsed, .convertToTime)
  rownames(mytime) <- c("STEP 0 (SVD decomposition)", "STEP 1 (variational algorithm)", "STEP 2 (summary statistics)", "STEP 3 (p0 estimation)", "STEP 4 (edge selection)", "overall")
  
  if(verbose){
    cat("\n\n")
    cat("Time (H:MM:SS):", .convertToTime(tps6[3]))
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
             time = mytime)
  
  return(out)
}
