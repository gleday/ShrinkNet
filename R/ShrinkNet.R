#' Gene Network Reconstruction using Global-Local Shrinkage Priors
#'
#' @param tX p by n data matrix of gene expression measurements
#' @param P matrix. Partition matrix for structured regularization.
#' @param globalShrink integer. Either 1 or 2. See Details.
#' @param p0 numeric. Prior estimate of the proportion of null hypotheses.
#' @param nsamp0 integer. If \code{p0}=\code{NULL}, number of randomly selected edges to estimate p0. See Details.
#' @param blfdr numeric. Bayesian analogue of the local false discovery rate used for edge selection. Value should be between 0 and 1. Default is 0.1.
#' @param maxNbEdges numeric. Maximum number of edges to select.
#' @param maxiter integer. Maximum number of iterations for the variational algorithm. Default is 100.
#' @param tol numeric. Represents the maximum relative convergence tolerance over the p variational lower bounds. Default is 0.001.
#' @param verbose logical. Should information on progress be printed?
#' @param standardize logical. Should the data be standardized? Default is TRUE.
#' @param sel logical. Whether edge selection should be carried out.
#' @param priorSig numeric. Vector of length 2 containg the shape and rate parameters of the gamma prior on the error precision.
#' @param weights matrix. Prior weights given to edges.
#' @param intercept logical. Whether an intercept should be included in the forward selection procedure.
#' @details
#' 
#' If \code{globalShrink}=1 then empirical Bayes for the global shrinkage prior is carried out using
#' fixed-point iterations as in Valpola and Honkela (2006). Otherwise, if \code{globalShrink}=2,
#' the approximate analytical solution of Leday et al (2016) is used.
#'
#' The user can specify using the argument \code{p0} a prior estimate of the proportion of edges not included in the graph,
#' which may have been determined from prior studies or public repositories. If \code{p0}=\code{NULL}, p0 is inferred from data.
#' In such case,  \code{nsamp0}=\code{NULL}, the proportion of null hypotheses p0 is estimated using Bayes factors calculated for all P=0.5*p*(p-1) edges (cf Leday et al., 2015).
#' When P is very large, it may me preferable to approximate p0 instead using a random subset of edges.
#' If \code{nsamp0} is an integer, then a random subset of size \code{nsamp0} is selected to estimate p0.
#' The default is \code{nsamp0}=\code{NULL} when p<=100 and \code{nsamp0}=1000 otherwise.
#' 
#' @return An object of class \code{\link{ShrinkNet-class}}
#' 
#' @author Gwenael G.R. Leday <gwenael.leday (at) mrc-bsu.cam.ac.uk>
#' 
#' @references Leday, G.G.R., de Gunst, M.C.M., Kpogbezan, G.B., van der Vaart, A.W., van Wieringen, W.N., and
#' van de Wiel, M.A. (2016). Gene network reconstruction using global-local shrinkage priors. To appear in The Annals of Applied Statistics.
#' @references Kpogbezan, G.B., van der Vaart, A.W., van Wieringen, W.N., Leday, G.G.R., and
#' van de Wiel, M.A. (2016). An empirical Bayes approach to network recovery using external knowledge. Submitted.
#' @references Valpola, H. and Honkela, A. (2006). Hyperparameter adaptation in variational Bayes for the gamma distribution.
#' Technical Report, Helsinki University of Technology.
#' 
#' @export
ShrinkNet <- function(tX, P=NULL , globalShrink=1, p0=NULL, nsamp0=NULL, blfdr=0.1, maxNbEdges=NULL, maxiter=100, tol=0.001, verbose=TRUE, standardize=TRUE, sel=TRUE, priorSig=c(0.001, 0.001), weights=NULL, intercept=TRUE){

  ##### Input checks
  if(is.matrix(tX)){
    if(any(is.na(tX))){
      stop("Missing values are not allowed")
    }
  }else{
    stop("tX is not a matrix")
  }
  if(!is.null(P)){
    if(is.matrix(P)){
      if(!isSymmetric(P)){
        stop("P is not symmetric")
      }
      if(nrow(P)!=nrow(tX)){
        stop("dimension of P must match the number of rows in tX")
      }
      diag(P) <- 0
      if(any(is.na(P))){
        stop("Missing values are not allowed")
      }
      Pvec <- P[upper.tri(P, diag=FALSE)]
      if(any(Pvec==0)){
        P <- P + 1
        diag(P) <- 0
        Pvec <- P[upper.tri(P, diag=FALSE)]
      }
      Pvals <- sort(unique(Pvec))
      K <- length(Pvals)
      if(K>2){
        stop("P cannot be partitioned in more than two groups")
      }
    }else{
      stop("P is not a matrix")
    }
  }
  if(!is.null(weights)){
    if(is.matrix(weights)){
      if(ncol(weights)==3){
        if(any(!as.vector(weights[,1:2])%in%c(1:nrow(tX)))){
          stop("wrong indexes in weights")
        }else{
          if(any(weights[,1]==weights[,2])){
            weights <- weights[!weights[,1]==weights[,2],]
          }else{
            if(any(weights[,3]<0) | any(weights[,3]>1)){
              stop("wrong edge weights")
            }else{
              tpw <- matrix(-1, nrow(tX), nrow(tX))
              matInd <- as.matrix(Matrix::sparseMatrix(i=weights[,1], j=weights[,2], x=TRUE, dims=c(nrow(tX), nrow(tX)), symmetric=TRUE))
              matVals <- as.matrix(Matrix::sparseMatrix(i=weights[,1], j=weights[,2], x=weights[,3], dims=c(nrow(tX), nrow(tX)), symmetric=TRUE))
              tpw[matInd] <- matVals[matInd]
              weights <- tpw
              rm(tpw, matInd, matVals)
            }
          }
        }
      }else{
        stop("weights must have three columns")
      }
    }else{
      stop("weights is not a matrix")
    }
  }else{
    weights <- matrix(-1, nrow(tX), nrow(tX))
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
  if(!is.logical(sel)){
    stop("sel is not a logical")
  }
  if(is.numeric(priorSig)){
    if(length(priorSig)!=2){
      stop("Two prior parameters exactly required in priorSig")
    }
    if(any(is.na(priorSig))){
      stop("priorSig contains NAs")
    }
    if(any(priorSig<=0)){
      stop("The two prior parameters in priorSig are not strictly positive")
    }
  }else{
    stop("priorSig is not a numeric")
  }
  tps <- proc.time()

  ##### Initialization
  aRand <- cSigma <- priorSig[1]
  bRand <- dSigma <- priorSig[2]
  
  ### If global shrinkage is not structured
  if(is.null(P)){
    
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
    eb <- HiddenVarAlgo(SVDs=allSVDs, tX=tX, aRand=aRand, bRand=bRand, cSigma=cSigma, dSigma=dSigma, maxiter=maxiter, globalShrink=globalShrink, tol=tol, verbose=verbose)
    
    tps2 <- proc.time() - tps - tps1
    
    ##### Calculate summary statistics from posteriors
    if(verbose) cat("STEP 2: Calculate summary statistics from posteriors... ")
    postSummaries <- sapply(1:nrow(tX), HiddenVarRidgeiGetKappa, SVDs=allSVDs, tX=tX, aRand=eb$parTau[nrow(eb$parTau),1], bRand=eb$parTau[nrow(eb$parTau),2], bRandStarInit=eb$allbRandStar, dSigmaStarInit=eb$alldSigmaStar, cSigma=cSigma, dSigma=dSigma, simplify=FALSE)
    matThres <- sapply(1:length(postSummaries), function(x){postSummaries[[x]][,1]}, simplify=TRUE)
    matBeta <- sapply(1:length(postSummaries), function(x){postSummaries[[x]][,2]}, simplify=TRUE)
    matThres <- (matThres + t(matThres))/2
    if(verbose) cat("DONE\n")
    
    tps3 <- proc.time() - tps - tps1 - tps2
    
    myP <- matrix(1, nrow(tX), nrow(tX))
    diag(myP) <- 0
    
    ##### Estimate p0
#     if(verbose) cat("STEP 3: Estimate p0... ")
#     if(sel){
#       if(is.null(p0)){
#         if(is.null(nsamp0)){
#           p0 <- HiddenEstimatep0(themat=matThres, tX=tX, cSigma=cSigma, dSigma=dSigma)
#           #myP <- matrix(1, nrow(tX), nrow(tX))
#           #diag(myP) <- 0
#           #resp0 <- HiddenEstimatep0P(themat=eb$matThres, tX=tX, P=myP, cSigma=cSigma, dSigma=dSigma)
#           #p0 <- resp0$p0
#         }else{
#           mat <- matThres
#           mat[upper.tri(mat)] <- 0
#           idx <- which(mat!=0, arr.ind=TRUE)
#           idx <- idx[sample(nrow(idx), nsamp0),]
#           allLogBFs <- t(apply(idx, 1, .edgeBFprime, themat=matThres, tX=tX, cSigma=cSigma, dSigma=dSigma))
#           #allLogBFs <- t(apply(idx-1, 1, HiddenEdgeBFprimeP, themat=eb$matThres, tX=tX, P=P, cSigma=cSigma, dSigma=dSigma))
#           p0 <- 1-mean(allLogBFs>0)
#         }
#       }else{
#         cat("(skipped)...")
#       }
#     }else{
#       p0 <- 0
#       cat("(skipped)...")
#     }
#     if(verbose) cat("DONE\n")
#     
#     tps4 <- proc.time() - tps - tps1 - tps2 - tps3
#     
#     ##### Edge selection using Bayesian local false discovery rate
#     if(verbose) cat("STEP 4: Edge selection... ")
#     if(sel){
#       resSel <- HiddenEdgeSelection(themat=matThres, tX=tX, p0=p0, lfdrcut=blfdr, maxNbEdges=maxNbEdges, cSigma=cSigma, dSigma=dSigma)
#       selGraph <- resSel$myGraph
#       logMaxBFs <- Matrix::Matrix(resSel$logMaxBFs, sparse=TRUE)
#       logMLSel <- resSel$sumLogML
#       nbedge <- sum(selGraph)/2
#       if(verbose){
#         cat("DONE\n\n")
#         cat("prior null probability p0 =", round(p0,5), "\n")
#         cat("", nbedge, " selected edges out of ", edgeTot, " (",round(100*nbedge/edgeTot, 2),"%)", " using blfdr = ", blfdr, sep="")
#       }
#     }else{
#       selGraph <- matrix(0, nrow(tX), nrow(tX))
#       logMaxBFs <- Matrix::Matrix(selGraph, sparse=TRUE)
#       logMLSel <- numeric(0)
#     }
#     tps5 <- proc.time() - tps - tps1 - tps2 - tps3 - tps4
#     tps6 <- proc.time() - tps
#     
#     ## Time 
#     mytime <- data.frame("elapsed"=c(tps1[3], tps2[3], tps3[3], tps4[3], tps5[3], tps6[3]))
#     mytime$"H:MM:SS" <- sapply(mytime$elapsed, .convertToTime)
#     rownames(mytime) <- c("STEP 0 (SVD decomposition)", "STEP 1 (variational algorithm)", "STEP 2 (summary statistics)", "STEP 3 (p0 estimation)", "STEP 4 (edge selection)", "overall")
#     
#     if(verbose){
#       cat("\n\n")
#       cat("Time (H:MM:SS):", .convertToTime(tps6[3]))
#       cat("\n\n")
#     }
#     
#     ## Output
#     myigraph <- igraph::graph.adjacency(selGraph, mode = "undirected")
#     if(!is.null(rownames(tX))){
#       myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=rownames(tX))
#     }
  }else{
    ##### Algo
    if(verbose) cat("STEP 1: Variational algorithm...\n")
    eb <- HiddenVarAlgoP(tX=tX, P=P, aRandInit=aRand, bRandInit=bRand, maxiter=maxiter, globalShrink=globalShrink, tol=tol, verbose=verbose, cSigma=cSigma, dSigma=dSigma)
    matThres <- eb$matThres
    matBeta <- eb$matBeta
    tps1 <- proc.time() - tps
    
    myP <- P
  }
  
  ##### Estimate p0
  if(is.null(P)){
    if(verbose) cat("STEP 3: Estimate p0... ")
  }else{
    if(verbose) cat("STEP 2: Estimate p0... ")
  }
  if(sel){
    if(is.null(p0)){
      if(is.null(nsamp0)){
        resp0 <- HiddenEstimatep0P(themat=matThres, tX=tX, P=myP, cSigma=cSigma, dSigma=dSigma, intercept=intercept)
        p0 <- resp0$p0
      }else{
        mat <- matThres
        mat[upper.tri(mat)] <- 0
        idx <- which(mat!=0, arr.ind=TRUE)
        idx <- idx[sample(nrow(idx), nsamp0),]
        allLogBFs <- t(apply(idx-1, 1, HiddenEdgeBFprimeP, themat=matThres, tX=tX, P=myP, cSigma=cSigma, dSigma=dSigma, intercept=intercept))
        p0 <- 1-mean(allLogBFs>0)
      }
    }else{
      if(verbose) cat("(skipped)...")
    }
  }else{
    p0 <- 0
    if(verbose) cat("(skipped)...")
  }
  if(verbose) cat("DONE\n")
  if(is.null(P)){
    tps4 <- proc.time() - tps - tps1 - tps2 - tps3
  }else{
    tps2 <- proc.time() - tps - tps1
  }
  
  ##### Edge selection using Bayesian local false discovery rate
  if(is.null(P)){
    if(verbose) cat("STEP 4: Edge selection... ")
  }else{
    if(verbose) cat("STEP 3: Edge selection... ")
  }
  if(sel){
    resSel <- HiddenEdgeSelectionP(themat=matThres, tX=tX, P=myP, weights=weights, p0=p0, lfdrcut=blfdr, maxNbEdges=maxNbEdges, cSigma=cSigma, dSigma=dSigma, intercept=intercept)
    selGraph <- resSel$myGraph
    logMaxBFs <- Matrix::Matrix(resSel$logMaxBFs, sparse=TRUE)
    logMLSel <- resSel$sumLogML
    nbedge <- sum(selGraph)/2
    if(verbose){
      cat("DONE\n\n")
      cat("prior null probability p0 =", round(p0,5), "\n")
      cat("", nbedge, " selected edges out of ", edgeTot, " (",round(100*nbedge/edgeTot, 2),"%)", " using blfdr = ", blfdr, sep="")
    }
  }else{
    selGraph <- matrix(0, nrow(tX), nrow(tX))
    logMaxBFs <- Matrix::Matrix(selGraph, sparse=TRUE)
    logMLSel <- numeric(0)
    cat("(skipped)... DONE\n\n")
  }
  if(is.null(P)){
    tps5 <- proc.time() - tps - tps1 - tps2 - tps3 - tps4
    tps6 <- proc.time() - tps
  }else{
    tps3 <- proc.time() - tps - tps1 - tps2
    tps4 <- proc.time() - tps
  }

  ## Time
  if(is.null(P)){
    mytime <- data.frame("elapsed"=c(tps1[3], tps2[3], tps3[3], tps4[3], tps5[3], tps6[3]))
    mytime$"H:MM:SS" <- sapply(mytime$elapsed, .convertToTime)
    rownames(mytime) <- c("STEP 0 (SVD decomposition)", "STEP 1 (variational algorithm)", "STEP 2 (summary statistics)", "STEP 3 (p0 estimation)", "STEP 4 (edge selection)", "overall")
    if(verbose){
      cat("\n\n")
      cat("Time (H:MM:SS):", .convertToTime(tps6[3]))
      cat("\n\n")
    }
  }else{
    mytime <- data.frame("elapsed"=c(tps1[3], tps2[3], tps3[3], tps4[3]))
    mytime$"H:MM:SS" <- sapply(mytime$elapsed, .convertToTime)
    rownames(mytime) <- c("STEP 1 (variational algorithm)", "STEP 2 (p0 estimation)", "STEP 3 (edge selection)", "overall")
    if(verbose){
      cat("\n\n")
      cat("Time (H:MM:SS):", .convertToTime(tps4[3]))
      cat("\n\n")
    }
  }
  
  ## Output
  myigraph <- igraph::graph.adjacency(selGraph, mode = "undirected")
  if(!is.null(rownames(tX))){
    myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=rownames(tX))
  }

  out <- new("ShrinkNet",
             graph = myigraph,
             kappa = matThres,
             beta = matBeta,
             p0 = p0,
             logMaxBFs = logMaxBFs,
             globalPrior = eb$parTau,
             allmargs = eb$allmargs,
             logMLSel = logMLSel,
             time = mytime)
  
  return(out)
}
