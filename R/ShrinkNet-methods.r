##########################
## Methods for "ShrinkNet"
##########################

#' @rdname ShrinkNet-class
#' @aliases print
#' @param x An object of class \code{ShrinkNet-class}
#' @param ... further arguments passed to or from other methods.
setMethod(
	f = "print",
	signature = "ShrinkNet",
	definition = function(x,...){
		cat("Class \"ShrinkNet\"\n")
		cat(nrow(x@kappa), " nodes", sep="")
		cat("\n")
	}
)

#' @rdname ShrinkNet-class
#' @aliases show
#' @param object An object of class \code{ShrinkNet-class}
setMethod(
	f = "show",
	signature = "ShrinkNet",
	definition = function(object) print(object)
)

#' @rdname ShrinkNet-class
#' @aliases summary
setMethod(
  f = "summary",
  signature = "ShrinkNet",
  definition = function(object,...){
    edgeTot <- 0.5*nrow(object@kappa)*(nrow(object@kappa)-1)
    cat("Class \"ShrinkNet\"\n\n")
    cat("UNDIRECTED GRAPH\n")
    cat("nodes:", vcount(object@graph), "\n")
    cat("edges:", ecount(object@graph), "\n")
    cat("density:", graph.density(object@graph), "\n\n")
    
    cat("Summary node degree:\n")
    print(summary(igraph::degree(object@graph)))
  }
)

#' @rdname ShrinkNet-class
#' @aliases adjacency
setMethod(
  f = "adjacency",
  signature = "ShrinkNet",
  definition = function(object){
    get.adjacency(object@graph)
  }
)

#' @rdname ShrinkNet-class
#' @aliases score
setMethod(
  f = "score",
  signature = "ShrinkNet",
  definition = function(object){
    object@kappa
  }
)

#' @rdname ShrinkNet-class
#' @aliases plotML
setMethod(
  f = "plotML",
  signature = "ShrinkNet",
  definition = function(object, ...){
    xx <- rowMeans(object@allmargs)
    plot(2:length(xx), xx[-1], xlab="iteration", ylab="average variational lower bound")
  }
)

#' @rdname ShrinkNet-class
#' @aliases plotPrior
setMethod(
  f = "plotPrior",
  signature = "ShrinkNet",
  definition = function(object, ...){
    xx <- rowMeans(object@allmargs)
    shapes <- object@globalPrior[,1]
    rates <- object@globalPrior[,2]
    pp <- length(shapes)
    x1 <- min(qgamma(0.001, shape=shapes[pp], rate=rates[pp]))
    x2 <- max(qgamma(0.999, shape=shapes[pp], rate=rates[pp]))
    xx <- seq(from=x1, to=x2, length=1000)
    plot(xx,dgamma(xx, shape=shapes[pp], rate=rates[pp]), xlim=c(x1,x2), type="l", xlab=expression(tau[j]^-2), ylab="probability", ...)
  }
)

#' @rdname ShrinkNet-class
#' @aliases plotGraph
setMethod(
  f = "plotGraph",
  signature = "ShrinkNet",
  definition = function(object, ...){
    plot(object@graph, ...)
  }
)
