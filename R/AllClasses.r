################
## All Classes
################

#' @name ShrinkNet-class
#' @aliases ShrinkNet-class
#'
#' @title Class ShrinkNet
#'
#' @description An S4 class representing the output of
#' the \code{\link{ShrinkNet}} function.
#'
#' @slot graph \code{\link[igraph]{igraph}} object
#' @slot kappa \code{matrix} of edge ranking scores
#' @slot beta \code{matrix} containing the posterior mean of
#' the SEM's regression coefficients (the ith column
#' corresponds to the ith regression equation).
#' @slot p0 \code{numeric}, estimate of the
#' proportion of null hypothesis.
#' @slot logMaxBFs \code{matrix} containing the log-max
#' Bayes factor for each edge (see Leday et. al (2015)).
#' @slot globalPrior \code{matrix} of the
#' (shape and rate) parameter values for the global shrinkage prior
#' (gamma distribution) at each iteration of the variational algorithm.
#' @slot allmargs \code{matrix} of variational lower bounds for
#' each regression equation in the SEM and each iteration of
#' the variational algorithm.
#' @slot time \code{numeric}, running time (in seconds).
#'
#' @method
#'  \item{print}{Print the object information}
#'  \item{show}{Print the object information}
#'
#' @importFrom "methods" "setClass" "representation"
#'
#' @author Gwenael G.R. Leday <gwenael.leday (at) mrc-bsu.cam.ac.uk>
#'
#' @references Leday, G.G.R., de Gunst, M.C.M., Kpogbezan, G.B.,
#' van der Vaart, A.W., van Wieringen, W.N., and
#' van de Wiel, M.A. (2015). Gene network reconstruction using
#' global-local shrinkage priors. The Annals of Applied Statistics.
#' 11 (2017), no. 1, 41--68.
#'
methods::setClass(
  "ShrinkNet",
  methods::representation(
    graph = "ANY",
    kappa = "matrix",
    beta = "matrix",
    p0 = "numeric",
    logMaxBFs = "ANY",
    globalPrior = "matrix",
    allmargs = "matrix",
    time = "data.frame"
  )
)
