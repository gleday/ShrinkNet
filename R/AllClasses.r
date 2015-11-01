################
## All Classes
################

#' @name ShrinkNet-class
#' @aliases ShrinkNet-class
#' 
#' @title Class ShrinkNet
#'
#' @description An S4 class representing the output of the \code{\link{ShrinkNet}} function.
#' 
#' @slot graph igraph. A graph object of class \code{\link[igraph]{igraph}}.
#' @slot kappa matrix. The matrix containing scores used to rank edges.
#' @slot p0 numeric. A number giving the estimate of the proportion of null hypothesis.
#' @slot globalPrior matrix. A matrix containing the (shape and rate) parameters of the global shrinkage prior
#' (gamma distribution) at each iteration of the variational algorithm.
#' @slot allmargs matrix. A matrix containing the values of the variational lower bounds for each regression
#' equation in the Bayesian SEM at each iteration of the variational algorithm.
#' @slot time numeric. Running time (in seconds) of ShrinkNet.
#' 
#' @method 
#'  \item{print}{Print the object information}
#'  \item{show}{Print the object information}
#'
#' @author Gwenael G.R. Leday <gwenael.leday (at) mrc-bsu.cam.ac.uk>
#' 
#' @references Leday, G.G.R., de Gunst, M.C.M., Kpogbezan, G.B., van der Vaart, A.W., van Wieringen, W.N., and
#' van de Wiel, M.A. (2015). Gene network reconstruction using global-local shrinkage priors. Submitted.
#'
setClass("ShrinkNet",
         representation(graph = "ANY",
                        kappa = "matrix",
                        p0 = "numeric",
                        globalPrior = "matrix",
                        allmargs = "matrix",
                        time = "numeric")
)
