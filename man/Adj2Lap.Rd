% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Adj2Lap.R
\name{Adj2Lap}
\alias{Adj2Lap}
\title{Compute graph Laplacian matrix from graph adjacency matrix}
\usage{
Adj2Lap(adj)
}
\arguments{
\item{adj}{graph adjacency matrix (squared symmetric matrix)}
}
\value{
graph Laplacian matrix
}
\description{
Compute graph Laplacian matrix from graph adjacency matrix
}
\examples{
# Define exemplary adjacency matrix
p1 <- 10
p2 <- 40
p <- p1 + p2
A <- matrix(rep(0, p * p), p, p)
A[1:p1, 1:p1] <- 1
A[(p1 + 1):p, (p1 + 1):p] <- 1
vizu.mat(A, "adjacency matrix")

# Compute corresponding Laplacian matrix
L <- Adj2Lap(A)
vizu.mat(L, "Laplacian matrix")
}
