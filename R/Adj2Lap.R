#' Compute graph Laplacian matrix from graph adjacency matrix 
#' 
#' @param adj graph adjacency matrix (squared symmetric matrix)
#' @return graph Laplacian matrix 
#' @examples 
#' # Define exemplary adjacency matrix
#' p1 <- 10
#' p2 <- 40
#' p <- p1 + p2
#' A <- matrix(rep(0, p * p), p, p)
#' A[1:p1, 1:p1] <- 1
#' A[(p1 + 1):p, (p1 + 1):p] <- 1
#' vizu.mat(A, "adjacency matrix")
#' 
#' # Compute corresponding Laplacian matrix
#' L <- Adj2Lap(A)
#' vizu.mat(L, "Laplacian matrix")
#' @export
#' 
Adj2Lap <- function(adj){
  D.diag <- apply(adj, MARGIN = 1, function(row) sum(row))
  D <- diag(D.diag, nrow = nrow(adj), ncol = ncol(adj))
  L <- D - adj
  return(L)
}