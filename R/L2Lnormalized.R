#' Compute normalized version of graph Laplacian matrix
#' 
#' @param L graph Laplcian matrix
#' @return normalized graph Laplacian matrix 
#' 
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
#' 
#' # Compute corresponding Laplacian matrix - normalized
#' L.norm <- L2L.normalized(L)
#' vizu.mat(L.norm, "L Laplacian matrix (normalized)")
#' 
#' @export
#' 
L2L.normalized <- function(L){
  L.d <- diag(L)
  L.normalized <- L
  p <- dim(L)[2]
  for (u in 1:p){
    for (v in 1:p){
      if (L.normalized[u, v] != 0){
        L.normalized[u, v] <- L.normalized[u,v] / sqrt(L.d[u] * L.d[v])
      }
    }
  }
  L.normalized.d <- rep(1, p)
  L.normalized.d[which(L.d == 0)] <- 0
  diag(L.normalized) <- L.normalized.d
  return(L.normalized)
}