################
## All Generics
################

setGeneric("adjacency", function(object) standardGeneric("adjacency"))
setGeneric("score", function(object) standardGeneric("score"))
setGeneric("plotML", function(object, ...) standardGeneric("plotML"))
setGeneric("plotPrior", function(object, ...) standardGeneric("plotPrior"))
setGeneric("plotGraph", function(object, ...) standardGeneric("plotGraph"))
setGeneric("topEdges", function(object, nb=10) standardGeneric("topEdges"))