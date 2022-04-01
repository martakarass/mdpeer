
#' Graph-constrained regression with addition of a small ridge term to handle the non-invertibility of a graph Laplacian matrix
#' 
#' @description 
#' Graph-constrained regression with addition of a diagonal matrix multiplied by a predefined (small) 
#' scalar to handle the non-invertibility of a graph Laplacian matrix (see: References). 
#' 
#' Bootstrap confidence intervals computation is available (not set as a default option). 
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}; typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling); 
#' **assumed to be already standarized** 
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling); 
#' **should contain colum of 1s if intercept is to be considered in a model** 
#' @param lambda.2 (small) scalar value of regularization parameter for diagonal matrix by adding which the \code{Q} matrix is corrected 
#' (note: correction is done *before* \eqn{\lambda_Q} regularization parameter value estimation; 
#' in other words: \eqn{\lambda_Q} estimation is done for the corrected \code{Q} matrix) 
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
#' \item{lambda.R}{\code{lambda.Q} * \code{lambda.2} value}
#' \item{lambda.2}{\code{lambda.2} supplied argument value}
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
#' riPEERc.out <- riPEERc(Q, y, Z, X)
#' plt.df <- data.frame(x = 1:p, y = riPEERc.out$b.est)
#' ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line() + labs("b estimates")
#' }
#' 
#' \dontrun{
#' # riPEERc with 0.95 bootstrap confidence intervals computation
#' riPEERc.out <- riPEERc(Q, y, Z, X, compute.boot.CI = TRUE, boot.R = 500)
#' plt.df <- data.frame(x = 1:p, y = riPEERc.out$b.est, 
#'                      lo = riPEERc.out$boot.CI[,1], 
#'                      up =  riPEERc.out$boot.CI[,2])
#' ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line() +  
#'   geom_ribbon(aes(ymin=lo, ymax=up), alpha = 0.3)
#' }
#' 
#' @references 
#' Karas, M., Brzyski, D., Dzemidzic, M., J., Kareken, D.A., Randolph, T.W., Harezlak, J. (2017). 
#' Brain connectivity-informed regularization methods for regression. doi: https://doi.org/10.1101/117945 
#' 
#' @import nlme
#' @import boot
#' @export
#' 
riPEERc <- function(Q, y, Z, X = NULL, 
                    lambda.2 = 0.001, 
                    compute.boot.CI = FALSE,
                    boot.R = 1000,
                    boot.conf = 0.95,
                    boot.set.seed = TRUE,
                    boot.parallel = "multicore",
                    boot.ncpus = 4,
                    verbose = TRUE){
  
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(y)[1] == dim(Z)[1])) stop("dim(y)[1] != dim(Z)[1]")
  if (!(dim(Q)[2] == dim(Z)[2])) stop("dim(Q)[2] != dim(Z)[2]")
  if (!is.null(X)) if(!(dim(y)[1] == dim(X)[1])) stop("dim(y)[1] != dim(X)[1]")
  
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  pcov <- ifelse(is.null(X), 0, dim(X)[2])
  
  # Regress out fixed effect coefficients (keep suffix ".wave" even if there is no fixed effect coefficients)
  if (pcov > 0){
    X.hash <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
    y.wave <- X.hash %*% y
    Z.wave <- X.hash %*% Z
  } else {
    y.wave <- y
    Z.wave <- Z
  }
  
  # Define LMM random effects matrix lme.Z
  Q_tilde <-  Q + lambda.2 * diag(p)
  Q_tilde.chol <- chol(Q_tilde)
  Q_tilde.chol.inv <- solve(Q_tilde.chol)
  lme.Z.wave <- Z.wave %*% Q_tilde.chol.inv
  
  # Fit LMM
  id.bd1 <- factor(rep(1, n))
  lme.fit <- lme(fixed = y.wave ~ 1, random = list(id.bd1 = pdIdent(~ lme.Z.wave - 1)), method = "REML") 
  
  # LMM-originated lambda parameter 
  sigma.eps <- lme.fit$sigma 
  sigma.u <- as.numeric(as.matrix(VarCorr(lme.fit))[1, 2])
  lambda.Q <- (sigma.eps^2)/(sigma.u^2)
  # Define lambda.R 
  lambda.R <- lambda.Q * lambda.2
  
  # Compute coefficient estimates
  b.est <- as.vector(solve(t(Z.wave) %*% Z.wave + lambda.Q * Q_tilde) %*% t(Z.wave) %*% y.wave)
  if (pcov > 0){
    beta.est <- as.vector(solve(t(X) %*% X) %*% t(X) %*% (y - Z %*% b.est))
  } else {
    beta.est <- c() 
  }
  
  # Compute bootstrap confidence intervals
  if (compute.boot.CI){
    if(verbose) message(paste0("riPEERc: computing bootstrap CI with ", boot.R, " replications..."))
    b.est.boot <- function(data, indices, lambda.Q, Q_tilde) {
      data.boot <- data[indices, ]
      y.wave.boot   <- as.matrix(data.boot[,1])
      Z.wave.boot  <- as.matrix(data.boot[,2:(ncol(data))])
      b.est.boot <- as.vector(solve(t(Z.wave.boot) %*% Z.wave.boot + lambda.Q * Q_tilde) %*% t(Z.wave.boot) %*% y.wave.boot)
      return(as.vector(b.est.boot))
    } 
    if (boot.set.seed) set.seed(1)
    boot.out <- boot(data = cbind(y.wave, Z.wave), statistic = b.est.boot,  R = boot.R, lambda.Q = lambda.Q, Q_tilde = Q_tilde, parallel = boot.parallel, ncpus = boot.ncpus)
    boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx){
      boot.ci.out <- (boot.ci(boot.out, type = "bca", index = idx, conf = boot.conf))$bca
      return(boot.ci.out[c(4,5)])
    })
    boot.CI <- do.call(rbind.data.frame, boot.CI.l)
    names(boot.CI) <- c("lower", "upper")
  } else {
    boot.CI <- NULL
  }
  
  res.list <- list(b.est    = b.est, 
                   beta.est = beta.est, 
                   lambda.Q = lambda.Q,
                   lambda.R = lambda.R,
                   lambda.2 = lambda.2,
                   boot.CI  = boot.CI)
  return(res.list)
}