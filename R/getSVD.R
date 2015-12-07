#' Convenience function for singular value decomposition
#'
#' @param ii integer. Gene index.
#' @param tX p by n matrix of gene expression.
#'
#' @details
#' The function returns the singular value decomposition of X_{ii}=UDV^T where
#' X_{ii} is the transpose of tX_{ii} which represents the matrix tX without the iith row.
#' 
#' @return A named list with the following elements:
#'  \item{u}{A matrix containing the left singular vectors.}
#'  \item{d}{A vector containing the singular values.}
#'  \item{v}{A matrix containing the right singular vectors.}
#'  
#' @export 
getSVD <- function(ii, tX) {
  .Call('ShrinkNet_getSVD', PACKAGE = 'ShrinkNet', ii, tX)
}
