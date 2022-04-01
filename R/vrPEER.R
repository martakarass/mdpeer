
#' Graph-constrained regression with variable-reduction procedure to handle the non-invertibility of 
#' a graph-originated penalty matrix
#' 
#' @description 
#' 
#' Graph-constrained regression with variable-reduction procedure to handle the non-invertibility of 
#' a graph-originated penalty matrix (see: References). 
#' 
#' Bootstrap confidence intervals computation is available (not set as a default option). 
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}; typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling); 
#' **assumed to be already standarized** 
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling); 
#' **should contain colum of 1s if intercept is to be considered in a model** 
#' @param sv.thr threshold value above which singular values of \code{Q} are considered "zeros" 
#' @param compute.boot.CI logical whether or not compute bootstrap confidence intervals for \eqn{b} regression coefficient estimates
#' @param boot.R number of bootstrap replications used in bootstrap confidence intervals computation
#' @param boot.conf confidence level assumed in bootstrap confidence intervals computation
#' @param boot.set.seed logical whether or not set seed in bootstrap confidence intervals computation
#' @param boot.parallel value of \code{parallel} argument in \code{boot} function in bootstrap confidence intervals computation
#' @param boot.ncpus value of \code{ncpus} argument in \code{boot} function in bootstrap confidence intervals computation
#' @param verbose logical whether or not set verbose mode (print out function execution messages)
#' @md 
#' 
#' @return 
#' \item{b.est}{vector of \eqn{b} coefficient estimates}
#' \item{beta.est}{vector of \eqn{\beta} coefficient estimates}
#' \item{lambda.Q}{\eqn{\lambda_Q} regularization parameter value}
#' \item{boot.CI}{data frame with two columns, \code{lower} and \code{upper}, containing, respectively, values of lower and upper bootstrap confidence intervals for \eqn{b} regression coefficient estimates}
#' 
#' @examples
#' set.seed(1234)
#' n <- 200
#' p1 <- 10
#' p2 <- 90
#' p <- p1 + p2
#' # Define graph adjacency matrix
#' A <- matrix(rep(0, p*p), nrow = p, ncol = p)
#' A[1:p1, 1:p1] <- 1
#' A[(p1+1):p, (p1+1):p] <- 1
#' L <- Adj2Lap(A)
#' # Define Q penalty matrix as graph Laplacian matrix normalized)
#' Q <- L2L.normalized(L)
#' # Define Z,X design matrices and aoutcome y
#' Z <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' b.true <- c(rep(1, p1), rep(0, p2))
#' X <- matrix(rnorm(n*3), nrow = n, ncol = 3)
#' beta.true <- runif(3)
#' intercept <- 0
#' eta <- intercept + Z %*% b.true + X %*% beta.true
#' R2 <- 0.5
#' sd.eps <- sqrt(var(eta) * (1 - R2) / R2)
#' error <- rnorm(n, sd = sd.eps)
#' y <- eta + error
#' 
#' \dontrun{
#' # run vrPEER 
#' vrPEER.out <- vrPEER(Q, y, Z, X)
#' plt.df <- data.frame(x = 1:p, 
#'                      y = vrPEER.out$b.est)
#' ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line()
#' }
#' 
#' \dontrun{
#' # run vrPEER with 0.95 confidence intrvals 
#' vrPEER.out <- vrPEER(Q, y, Z, X, compute.boot.CI = TRUE, boot.R = 500)
#' plt.df <- data.frame(x = 1:p, 
#'                      y = vrPEER.out$b.est, 
#'                      lo = vrPEER.out$boot.CI[,1], 
#'                      up =  vrPEER.out$boot.CI[,2])
#' ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line() +  
#'   geom_ribbon(aes(ymin=lo, ymax=up), alpha = 0.3)
#' }
#' 
#' @references 
#' Karas, M., Brzyski, D., Dzemidzic, M., J., Kareken, D.A., Randolph, T.W., Harezlak, J. (2017). 
#' Brain connectivity-informed regularization methods for regression. doi: https://doi.org/10.1101/117945 
#' 
#' @import boot
#' @import nlme
#' @export
#' 
vrPEER <- function(Q, y, Z, X = NULL, 
                   sv.thr = 1e-5, 
                   compute.boot.CI = FALSE,
                   boot.R = 1000,
                   boot.conf = 0.95,
                   boot.set.seed = TRUE,
                   boot.parallel = "multicore",
                   boot.ncpus = 4,
                   verbose = TRUE){

  # Check for data objects dimensions consistency 
  if (dim(y)[1] != dim(Z)[1]) stop("dim(y)[1] != dim(Z)[1]")
  if (dim(Q)[2] != dim(Z)[2]) stop("dim(L)[2] != dim(Z)[2]")
  if (!is.null(X)) if(dim(y)[1] != dim(X)[1]) stop("dim(y)[1] != dim(X)[1]")

  # Construct data objects used in variable-reduction approach
  pcov <- ifelse(is.null(X), 0, dim(X)[2])
  Z.ext <- cbind(X, Z) 
  Q.ext <- as.matrix(magic::adiag(matrix(0, nrow = pcov, ncol = pcov), Q))  # Q extension with 0-filled matrix to reflect penalty for predictors we do not want to penalize in estimation
  n <- dim(Z.ext)[1]
  p <- dim(Z.ext)[2]
  
  # SVD of Q 
  Q.ext.svd <- svd(Q.ext)
  U <- Q.ext.svd$u
  Sigma <- diag(Q.ext.svd$d)
  r <- length(which(diag(Sigma) < sv.thr))  # how many singular values of Sigma are "zeros"
  k <- p - r # how many singular values of Sigma are non-"zeros" (are OK)
  
  if (r == 0){
    # Define LMM random effects matrix lme.Z
    Q.chol <- chol(Q)
    Q.chol.inv <- solve(Q.chol)
    lme.Z <- Z %*% Q.chol.inv
    # Fit LMM
    id.bd1 <- factor(rep(1, n))
    lme.fit <- lme(fixed = y ~ 1, random = list(id.bd1 = pdIdent(~ lme.Z - 1)), method = "REML") 
  } else {
    # Procedure step (1): transform PEER to weighted Ridge
    Sigma.wave <- diag(diag(Sigma)[1:k])  # reduce Sigma to keep only "non-zeros" singular values 
    Z.wave <- Z.ext %*% U
    X.wave <- as.matrix(Z.wave[, (k+1):p])  # fixed effects matrix 
    Z.wave.red <- Z.wave[, 1:k]  # reduce Z.wave to keep only "non-zeros" singular values 
    # Procedure step (2): transform weighted Ridge to Ridge 
    b <- diag(sqrt(diag(Sigma.wave)))
    b.inv <- solve(b)
    Z.wave.red.wave <- Z.wave.red %*% b.inv
    # Fit LMM model
    id.bd1 <- factor(rep(1, n))
    lme.fit <- lme(fixed = y ~ X.wave + 1, random = list(id.bd1 = pdIdent(~ Z.wave.red.wave - 1)), method = "REML") 
  }
  
  # LMM-originated lambda parameter
  sigma.eps <- lme.fit$sigma 
  sigma.u <- as.numeric(as.matrix(VarCorr(lme.fit))[1, 2])
  lambda.Q <- (sigma.eps^2)/(sigma.u^2)
  
  # Compute coefficient estimates
  beta.b.est <- as.vector(solve(t(Z.ext) %*% Z.ext + lambda.Q * Q.ext) %*% t(Z.ext) %*% y)
  if (is.null(colnames(Z.ext))) names(beta.b.est) <- colnames(Z.ext)
  if (pcov > 0){
    b.est <- beta.b.est[(pcov+1):(pcov+dim(Z)[2])]
    beta.est <- beta.b.est[1:pcov]
  } else {
    b.est <- beta.b.est
    beta.est <- c()
  }
  
  # Compute bootstrap confidence intervals
  if (compute.boot.CI){
    if(verbose) message(paste0("vrPEER: computing bootstrap CI with ", boot.R, " replications..."))
    b.est.boot <- function(data, indices, lambda.Q, Q.ext, pX, pZ) {
      data.boot <- data[indices, ]
      y.boot      <- as.matrix(data.boot[,1])
      Z.ext.boot  <- as.matrix(data.boot[,2:(ncol(data))])
      beta.b.est.boot <- as.vector(solve(t(Z.ext.boot) %*% Z.ext.boot + lambda.Q * Q.ext) %*% t(Z.ext.boot) %*% y.boot)
      if (pX > 0){
        b.est.boot <- beta.b.est.boot[(pX+1):(pX+pZ)]
      } else {
        b.est.boot <- beta.b.est.boot
      }
      return(as.vector(b.est.boot))
    } 
    if (boot.set.seed) set.seed(1)
    boot.out <- boot(data = cbind(y, Z.ext), statistic = b.est.boot, 
                     R = boot.R, lambda.Q = lambda.Q, Q.ext = Q.ext, pX = pcov, pZ = dim(Z)[2], 
                     parallel = boot.parallel, ncpus = boot.ncpus)
    boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx){
      boot.ci.out <- (boot.ci(boot.out, type = "bca", index = idx, conf = boot.conf))$bca
      return(boot.ci.out[c(4,5)])
    })
    boot.CI <- do.call(rbind.data.frame, boot.CI.l)
    names(boot.CI) <- c("lower", "upper")
  } else {
    boot.CI <- NULL
  }
  
  res.list <- list(b.est = b.est, 
                   beta.est = beta.est, 
                   lambda.Q = lambda.Q,
                   boot.CI  = boot.CI)
  return(res.list)
}