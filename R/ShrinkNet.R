#' Gene network reconstruction using global-local shrinkage priors
#'
#' @param tX \code{matrix}. Gene expression data with \eqn{p} rows
#' (markers) and \eqn{n} columns (obervations).
#' @param globalShrink \code{integer}. Type of global shrinkage,
#' equal to 1 or 2. See Details.
#' @param nsamp0 \code{integer}. Number of randomly selected
#' edges used to estimate \eqn{p_0}. See Details.
#' @param blfdr \code{numeric}. Bayesian analogue of the
#' local false discovery rate used for
#' edge selection (value between 0 and 1). Default is 0.1.
#' @param maxNbEdges \code{numeric}. Maximum number of edges to select.
#' @param maxiter \code{integer}. Maximum number of iterations
#' for the variational algorithm. Default is 100.
#' @param tol \code{numeric}. Maximum relative convergence tolerance
#' over the \eqn{p} variational lower bounds. Default is 0.001.
#' @param verbose \code{logical}. Should information on progress be printed?
#' @param standardize \code{logical}. Should the data be
#' standardized? Default is TRUE.
#' @details
#'
#' When \code{globalShrink} = 1, empirical Bayes for the global
#' shrinkage prior is carried out using fixed-point iterations as
#' in Valpola and Honkela (2006). Otherwise, when \code{globalShrink} = 2,
#' the approximate analytical solution of Leday et al (2015) is used.
#'
#' When \code{nsamp0} = \code{NULL} (default), the proportion of null
#' hypotheses \eqn{p_0} is estimated using the Bayes factors
#' calculated for all \eqn{P = 0.5 * p * (p - 1)} edges,
#' as in Leday et al (2015). However, when \eqn{P} is very large,
#' it may be preferable (for computational reasons) to
#' estimate \eqn{p_0} from a random subset of edges.
#' The argument \code{nsamp0} allows to specify the size of
#' the random subset. Currently, the default is
#' \code{nsamp0} = \code{NULL} when \eqn{p <= 100} and
#' \code{nsamp0} = 1000 otherwise.
#'
#' @useDynLib ShrinkNet
#' @import Rcpp
#' @importFrom "Matrix" "Matrix"
#' @importFrom "igraph" "set_vertex_attr" "graph.adjacency"
#' @importFrom "methods" "new"
#' @importFrom "assertthat" "assert_that" "not_empty" "noNA" "is.count"
#' @importFrom "dplyr" "between"
#'
#' @return An object of class \code{\link{ShrinkNet-class}}
#'
#' @author Gwenael G.R. Leday
#'
#' @references Leday, G.G.R., de Gunst, M.C.M., Kpogbezan, G.B.,
#' van der Vaart, A.W., van Wieringen, W.N., and
#' van de Wiel, M.A. (2015). Gene network reconstruction using
#' global-local shrinkage priors. The Annals of Applied Statistics.
#' 11 (2017), no. 1, 41--68.
#'
#' @export
ShrinkNet <- function(
  tX,
  globalShrink = 1,
  nsamp0 = NULL,
  blfdr = 0.1,
  maxNbEdges = NULL,
  maxiter = 100,
  tol = 0.001,
  verbose = TRUE,
  standardize = TRUE
) {

  # check tX
  assertthat::assert_that(
    is.matrix(tX),
    assertthat::not_empty(tX),
    assertthat::noNA(tX),
    all(is.finite(tX))
  )

  # check globalShrink
  assertthat::assert_that(
    assertthat::is.count(globalShrink),
    is.finite(globalShrink),
    globalShrink %in% c(1, 2)
  )

  # number of edges
  edgeTot <- 0.5 * nrow(tX) * (nrow(tX) - 1)

  # check nsamp0
  nsamp0 <- ifelse(
    is.null(nsamp0),
    ifelse(nrow(tX) > 100, 1000, edgeTot),
    nsamp0
  )
  assertthat::assert_that(
    assertthat::is.count(nsamp0),
    is.finite(nsamp0)
  )
  nsamp0 <- pmin(nsamp0, edgeTot)
  if (nsamp0 < 1000) {
    warning("nsamp0 (< 1000) may be too low to
            obtain a reasonable estimate of p0")
  }

  # maxiter
  assertthat::assert_that(
    assertthat::is.count(maxiter),
    is.finite(maxiter)
  )

  # blfdr
  assertthat::assert_that(
    assertthat::is.number(blfdr),
    is.finite(blfdr),
    dplyr::between(blfdr, 0, edgeTot),
    msg = "blfdr must be a number between 0 and 1"
  )

  # maxNbEdges
  maxNbEdges <- ifelse(is.null(maxNbEdges), 0, maxNbEdges)
  assertthat::assert_that(
    assertthat::is.count(maxNbEdges),
    is.finite(maxNbEdges),
    dplyr::between(maxNbEdges, 1, edgeTot),
    msg = paste0(
      "maxNbEdges must be a number between 1 and ",
      edgeTot
    )
  )

  # verbose
  assertthat::assert_that(
    is.logical(verbose),
    assertthat::is.scalar(verbose),
    assertthat::not_empty(verbose)
  )

  # standardize
  assertthat::assert_that(
    is.logical(standardize),
    assertthat::is.scalar(standardize),
    assertthat::not_empty(standardize)
  )
  if (standardize) {
    tX <- t(scale(t(tX), center = TRUE, scale = TRUE))
  } else {
    warning("Input data have not been standardized")
  }

  tps <- proc.time()

  ##### Initialization
  aRand <- 0.001
  bRand <- 0.001

  ##### Data preparation
  if (verbose) {
    cat("\n")
    cat("STEP 0: SVD computations... ")
  }

  allSVDs <- sapply(seq_len(nrow(tX)), getSVD, tX = tX, simplify = FALSE)

  if (verbose) cat("DONE\n")
  tps1 <- proc.time() - tps

  ##### Algo
  if (verbose) cat("STEP 1: Variational algorithm...\n")

  eb <- HiddenVarAlgo(
    SVDs = allSVDs,
    tX = tX,
    aRand = aRand,
    bRand = bRand,
    maxiter = maxiter,
    globalShrink = globalShrink,
    tol = tol,
    verbose = verbose
  )

  tps2 <- proc.time() - tps - tps1

  ##### Calculate summary statistics from posteriors
  if (verbose) cat("STEP 2: Calculate summary statistics from posteriors... ")

  postSummaries <- sapply(
    seq_len(nrow(tX)),
    HiddenVarRidgeiGetKappa,
    SVDs = allSVDs,
    tX = tX,
    aRand = eb$parTau[nrow(eb$parTau), 1],
    bRand = eb$parTau[nrow(eb$parTau), 2],
    bRandStarInit = eb$allbRandStar,
    dSigmaStarInit = eb$alldSigmaStar,
    simplify = FALSE
  )
  matThres <- sapply(
    seq_len(length(postSummaries)),
    function(x) {
      postSummaries[[x]][, 1]
    },
    simplify = TRUE
  )
  matBeta <- sapply(
    seq_len(length(postSummaries)),
    function(x) {
      postSummaries[[x]][, 2]
    },
    simplify = TRUE
  )
  matThres <- (matThres + t(matThres)) / 2

  if (verbose) cat("DONE\n")
  tps3 <- proc.time() - tps - tps1 - tps2

  ##### Estimate p0
  if (verbose) cat("STEP 3: Estimate p0... ")

  if (is.null(nsamp0)) {
    p0 <- HiddenEstimatep0(themat = matThres, tX = tX)
  } else {
    mat <- matThres
    mat[upper.tri(mat)] <- 0
    idx <- which(mat != 0, arr.ind = TRUE)
    idx <- idx[sample(nrow(idx), nsamp0), ]
    allLogBFs <- t(apply(
      idx,
      1,
      .edgeBFprime,
      themat = matThres,
      tX = tX
    ))
    p0 <- 1 - mean(allLogBFs > 0)
  }

  if (verbose) cat("DONE\n")
  tps4 <- proc.time() - tps - tps1 - tps2 - tps3

  ##### Edge selection using Bayesian local false discovery rate
  if (verbose) cat("STEP 4: Edge selection... ")

  resSel <- HiddenEdgeSelection(
    themat = matThres,
    tX = tX,
    p0 = p0,
    lfdrcut = blfdr,
    maxNbEdges = maxNbEdges
  )
  selGraph <- resSel$myGraph
  logMaxBFs <- Matrix::Matrix(resSel$logMaxBFs, sparse = TRUE)
  nbedge <- sum(selGraph) / 2

  if (verbose) {
    cat("DONE\n\n")
    cat("prior null probability p0 =", round(p0, 5), "\n")
    cat(
      "",
      nbedge,
      " selected edges out of ",
      edgeTot,
      " (", round(100 * nbedge / edgeTot, 2), "%)",
      " using blfdr = ",
      blfdr,
      sep = ""
    )
  }
  tps5 <- proc.time() - tps - tps1 - tps2 - tps3 - tps4
  tps6 <- proc.time() - tps

  ## Time
  mytime <- data.frame(
    "elapsed" = c(
      tps1[3], tps2[3], tps3[3], tps4[3], tps5[3], tps6[3]
    )
  )
  mytime$"H:MM:SS" <- sapply(mytime$elapsed, .convert2time)
  rownames(mytime) <- c(
    "STEP 0 (SVD decomposition)",
    "STEP 1 (variational algorithm)",
    "STEP 2 (summary statistics)",
    "STEP 3 (p0 estimation)",
    "STEP 4 (edge selection)",
    "overall"
  )

  if (verbose) {
    cat("\n\n")
    cat("Time (H:MM:SS):", .convert2time(tps6[3]))
    cat("\n\n")
  }

  ## Output
  myigraph <- igraph::graph.adjacency(selGraph, mode = "undirected")
  if (!is.null(rownames(tX))) {
    myigraph <- igraph::set_vertex_attr(myigraph, "name", value = rownames(tX))
  }

  methods::new(
    "ShrinkNet",
    graph = myigraph,
    kappa = matThres,
    beta = matBeta,
    p0 = p0,
    logMaxBFs = logMaxBFs,
    globalPrior = eb$parTau,
    allmargs = eb$allmargs,
    time = mytime
  )
}
