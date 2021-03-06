% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vrPEER.R
\name{vrPEER}
\alias{vrPEER}
\title{Graph-constrained regression with variable-reduction procedure to handle the non-invertibility of
a graph-originated penalty matrix}
\usage{
vrPEER(Q, y, Z, X = NULL, sv.thr = 1e-05, compute.boot.CI = FALSE,
  boot.R = 1000, boot.conf = 0.95, boot.set.seed = TRUE,
  boot.parallel = "multicore", boot.ncpus = 4, verbose = TRUE)
}
\arguments{
\item{Q}{graph-originated penalty matrix \eqn{(p \times p)}; typically: a graph Laplacian matrix}

\item{y}{response values matrix \eqn{(n \times 1)}}

\item{Z}{design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling);
\strong{assumed to be already standarized}}

\item{X}{design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling);
\strong{should contain colum of 1s if intercept is to be considered in a model}}

\item{sv.thr}{threshold value above which singular values of \code{Q} are considered "zeros"}

\item{compute.boot.CI}{logical whether or not compute bootstrap confidence intervals for \eqn{b} regression coefficient estimates}

\item{boot.R}{number of bootstrap replications used in bootstrap confidence intervals computation}

\item{boot.conf}{confidence level assumed in bootstrap confidence intervals computation}

\item{boot.set.seed}{logical whether or not set seed in bootstrap confidence intervals computation}

\item{boot.parallel}{value of \code{parallel} argument in \code{boot} function in bootstrap confidence intervals computation}

\item{boot.ncpus}{value of \code{ncpus} argument in \code{boot} function in bootstrap confidence intervals computation}

\item{verbose}{logical whether or not set verbose mode (print out function execution messages)}
}
\value{
\item{b.est}{vector of \eqn{b} coefficient estimates}
\item{beta.est}{vector of \eqn{\beta} coefficient estimates}
\item{lambda.Q}{\eqn{\lambda_Q} regularization parameter value}
\item{boot.CI}{data frame with two columns, \code{lower} and \code{upper}, containing, respectively, values of lower and upper bootstrap confidence intervals for \eqn{b} regression coefficient estimates}
}
\description{
Graph-constrained regression with variable-reduction procedure to handle the non-invertibility of
a graph-originated penalty matrix (see: References).

Bootstrap confidence intervals computation is available (not set as a default option).
}
\examples{
set.seed(1234)
n <- 200
p1 <- 10
p2 <- 90
p <- p1 + p2
# Define graph adjacency matrix
A <- matrix(rep(0, p*p), nrow = p, ncol = p)
A[1:p1, 1:p1] <- 1
A[(p1+1):p, (p1+1):p] <- 1
L <- Adj2Lap(A)
# Define Q penalty matrix as graph Laplacian matrix normalized)
Q <- L2L.normalized(L)
# Define Z,X design matrices and aoutcome y
Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
b.true <- c(rep(1, p1), rep(0, p2))
X <- matrix(rnorm(n*3), nrow = n, ncol = 3)
beta.true <- runif(3)
intercept <- 0
eta <- intercept + Z \%*\% b.true + X \%*\% beta.true
R2 <- 0.5
sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
error <- rnorm(n, sd = sd.eps)
y <- eta + error

\dontrun{
# run vrPEER 
vrPEER.out <- vrPEER(Q, y, Z, X)
plt.df <- data.frame(x = 1:p, 
                     y = vrPEER.out$b.est)
ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line()
}

\dontrun{
# run vrPEER with 0.95 confidence intrvals 
vrPEER.out <- vrPEER(Q, y, Z, X, compute.boot.CI = TRUE, boot.R = 500)
plt.df <- data.frame(x = 1:p, 
                     y = vrPEER.out$b.est, 
                     lo = vrPEER.out$boot.CI[,1], 
                     up =  vrPEER.out$boot.CI[,2])
ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line() +  
  geom_ribbon(aes(ymin=lo, ymax=up), alpha = 0.3)
}

}
\references{
Karas, M., Brzyski, D., Dzemidzic, M., J., Kareken, D.A., Randolph, T.W., Harezlak, J. (2017).
Brain connectivity-informed regularization methods for regression. doi: https://doi.org/10.1101/117945
}
