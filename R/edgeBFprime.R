# Internal function
# Wrapper for ShrinkNet_HiddenEdgeBFprime

# Author: Gwenael G.R. Leday

.edgeBFprime <- function(idx, themat, tX){
  .Call('ShrinkNet_HiddenEdgeBFprime', PACKAGE = 'ShrinkNet', idx, themat, tX)
}