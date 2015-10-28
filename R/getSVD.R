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
#'  \item{myF}{A matrix representing the crossproduct u*diag(D).}
#'  \item{FTF}{A matrix representing the crossproduct t(myF)*myF.}
#'  \item{myy}{A vector containing expression values for gene ii.}
#'  \item{myX}{A n by p-1 matrix containing expression values for all genes other than ii.}
#'  
#' @export 
getSVD <- function(ii, tX) {
  .Call('ShrinkNet_getSVD', PACKAGE = 'ShrinkNet', ii, tX)
}
