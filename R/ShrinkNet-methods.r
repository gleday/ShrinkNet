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
    cat("nodes:", igraph::vcount(object@graph), "\n")
    cat("edges:", igraph::ecount(object@graph), "\n")
    cat("density:", igraph::graph.density(object@graph), "\n\n")
    
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
    igraph::get.adjacency(object@graph)
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

#' @rdname ShrinkNet-class
#' @aliases listEdges
setMethod(
  f = "listEdges",
  signature = "ShrinkNet",
  definition = function(object){
    mat <- igraph::get.edgelist(object@graph, names = TRUE)
    mat <- as.data.frame(mat, stringsAsFactors=FALSE)
    labs <- igraph::get.vertex.attribute(object@graph, "name")
    indx <- t(apply(mat, 1, function(x){c(which(labs==x[1]), which(labs==x[2]))}))
    mat <- cbind(mat, indx)
    colnames(mat) <- c("node1","node2","index1","index2")
    mat$score <- apply(indx, 1, function(x){object@kappa[x[1],x[2]]})
    mat <- mat[order(mat$score, decreasing=TRUE),]
    rownames(mat) <- NULL
    mat
  }
)

#' @rdname ShrinkNet-class
#' @aliases topEdges
#' @param nb maximum number of edges (or nodes, depending on the function it applies to) to return, Default is 20.
setMethod(
  f = "topEdges",
  signature = "ShrinkNet",
  definition = function(object, nb=20){
    themat <- score(object)
    themat[lower.tri(themat, diag=TRUE)] <- NA
    indx <- which(!is.na(themat), arr.ind = TRUE)
    myscores <- apply(indx, 1, function(x){object@kappa[x[1],x[2]]})
    labs <- igraph::get.vertex.attribute(object@graph, "name")
    nodes <- t(apply(indx, 1, function(x){c(labs[x[1]], labs[x[2]])}))
    mat <- data.frame(nodes, indx, myscores, stringsAsFactors=FALSE)
    colnames(mat) <- c("node1","node2","index1","index2","score")
    mat <- mat[order(mat$score, decreasing=TRUE),]
    rownames(mat) <- NULL
    mat[1:nb,]
  }
)

#' @rdname ShrinkNet-class
#' @aliases topDegree
setMethod(
  f = "topDegree",
  signature = "ShrinkNet",
  definition = function(object, nb=20){
    vals <- igraph::degree(object@graph)
    vals <- vals[order(vals, decreasing=TRUE)]
    mat <- matrix(vals[1:nb])
    rownames(mat) <- names(vals)[1:nb]
    colnames(mat) <- "degree"
    mat
  }
)


