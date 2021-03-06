#' mdpeer: Methods for graph-constrained regression with enhanced regularization parameters selection
#'
#' Provides graph-constrained
#' regression methods in which regularization parameters 
#' are selected automatically via estimation of equivalent Linear
#' Mixed Model formulation.  
#' 'riPEER' (ridgified Partially Empirical Eigenvectors for
#' Regression) method employs a penalty term being
#' a linear combination of graph-originated and ridge-originated penalty terms,
#' whose two regularization parameters are ML estimators from corresponding Linear
#' Mixed Model solution; a graph-originated penalty
#' term allows imposing similarity between coefficients based on graph information
#' given whereas additional ridge-originated penalty term facilitates parameters
#' estimation: it reduces computational issues arising from singularity in a graph-originated 
#' penalty matrix and yields plausible results in situations when graph
#' information is not informative. 'riPEERc' (ridgified Partially Empirical Eigenvectors for
#' Regression with constant) method utilizes addition of a diagonal matrix multiplied by a predefined (small) 
#' scalar to handle the non-invertibility of a graph Laplacian matrix. 'vrPEER' 
#' (variable reducted PEER) method performs variable-reduction procedure to handle 
#' the non-invertibility of a graph Laplacian matrix. 
#'
#' @docType package
#' @name mdpeer
NULL