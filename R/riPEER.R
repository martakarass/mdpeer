

# Wrapper to apply default values of selected arguments of some other function
# 
# @param func function to be wrapped
# 
# @return wrapped function
# 
applyDefaults <- function(func, ...){
  function(x) func(x, ...)
}




# Compute riPEER objective function value
#
# @param lambda function parameter
# @param S function parameter
# @param y function parameter
# @param Z function parameter
# 
# @return objective function value
# 
riPEER.obj.fun <- function(lambda, S, y, Z){
  lambda.Q <- lambda[1]
  lambda.R <- lambda[2]
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  D <- lambda.Q * S + lambda.R * diag(p)
  val <- n * log(sum(y^2) -  t(y) %*% Z %*% solve(D + t(Z) %*% Z) %*% t(Z) %*% y) + log(det((D + t(Z) %*% Z) %*% solve(D)))
  return(val)
}




# Define riPEER rootSolve algorithm stationary conditions function(s) 
#
# @param lambda function parameter
# @param S function parameter
# @param y function parameter
# @param Z function parameter
# 
# @return vector with stationary conditions function(s) values
# 
rootSolve.riPEER <- function(x, S, y, Z) {
  lambda <- abs(x)
  lambda.Q <- lambda[1]
  lambda.R <- lambda[2]
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  v <- t(Z) %*% y
  D_lambda <- lambda.Q * S + lambda.R * diag(p)
  W <- solve(D_lambda + t(Z) %*% Z)
  K <- 1/(sum(y^2) - t(v) %*% W %*% v)
  gr1 <- n * t(v) %*% W %*% S %*% W %*% v %*% K + psych::tr(W %*% S - solve(D_lambda) %*% S)
  gr2 <- n * t(v) %*% W %*% W %*% v %*% K + psych::tr(W - solve(D_lambda))
  c(F1 = gr1, F2 = gr2)
}




# Define riPEER rootSolve algorithm stationary conditions function(s) - 
# optimizing lambda.R only (lambda.Q=0 fixed)
#
# @param lambda function parameter
# @param S function parameter
# @param y function parameter
# @param Z function parameter
# 
# @return vector with stationary conditions function(s) values
#
rootSolve.Q0.riPEER <- function(x, S, y, Z) {
  lambda.Q <- 0
  lambda.R <- abs(x)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  v <- t(Z) %*% y
  D_lambda <- lambda.Q * S + lambda.R * diag(p)
  W <- solve(D_lambda + t(Z) %*% Z)
  K <- 1/(sum(y^2) - t(v) %*% W %*% v)
  gr2 <- n * t(v) %*% W %*% W %*% v %*% K + psych::tr(W - solve(D_lambda))
  c(F2 = gr2)
}




# Optimize lambda=(lambda.Q,lambda.R) regularization parameter values.
#
# @param S function parameter
# @param y function parameter
# @param Z function parameter
# 
# @return vector with stationary conditions function(s) values
#
riPEER.opt.lambda <- function(S, y, Z, 
                              optim.metod, 
                              sbplx.x0, sbplx.lambda.lo, sbplx.lambda.up,
                              rootSolve.x0, rootSolve.Q0.x0,
                              verbose) {
  obj.fn <- applyDefaults(riPEER.obj.fun, S = S, y = y, Z = Z)
  out <- tryCatch({
    if (optim.metod == "sbplx"){
      sbplx.out             <- sbplx(x0 = sbplx.x0, fn = obj.fn, lower = sbplx.lambda.lo, upper = sbplx.lambda.up) 
      lambda                <- sbplx.out$par
      obj.fn.val            <- obj.fn(lambda)
      list(lambda = lambda, obj.fn.val = obj.fn.val)
    } else if (optim.metod == "rootSolve") {
      multiroot.fn          <- applyDefaults(rootSolve.riPEER, S = S, y = y, Z = Z)
      multiroot.Q0.fn       <- applyDefaults(rootSolve.Q0.riPEER, S = S, y = y, Z = Z)
      # optimize lambdaQ, lambdaR
      multiroot.out.tmp     <- multiroot(f = multiroot.fn, start = rootSolve.x0)
      lambda.riPEER.min     <- abs(multiroot.out.tmp$root)
      obj.fun.riPPER.min    <- obj.fn(lambda.riPEER.min)
      # if (verbose) message(paste0("riPEER:   obj.fn: ", obj.fun.riPPER.min, ",   lambda: ", paste0(lambda.riPEER.min, collapse = ", ")))
      # lambdaQ=0, optimize lambdaR
      multiroot.Q0.out.tmp  <- multiroot(f = multiroot.Q0.fn, start = rootSolve.Q0.x0)
      lambda.riPEER.Q0.min  <- c(0, abs(multiroot.Q0.out.tmp$root))
      obj.fun.riPPER.Q0.min <- obj.fn(lambda.riPEER.Q0.min)
      # if (verbose) message(paste0("riPEER.Q0: obj.fn: ", obj.fun.riPPER.Q0.min, ",   lambda: ", paste0(lambda.riPEER.Q0.min, collapse = ", ")))
      # select minimum between riPEER, riPEER.Q0
      if (obj.fun.riPPER.min < obj.fun.riPPER.Q0.min){
        lambda.OPT  <- lambda.riPEER.min
        obj.fun.OPT <- obj.fun.riPPER.min
      } else {
        lambda.OPT  <- lambda.riPEER.Q0.min
        obj.fun.OPT <- obj.fun.riPPER.Q0.min
      }
      return(list(lambda = lambda.OPT, obj.fn.val = obj.fun.OPT)) 
    } else {
      NULL
    }
  }, error = function(e) {
    message(e)
    NULL
  })
  return(out)
}




# ------------------------------------------------------------------------------

#' Graph-constrained regression with penalty term being a linear combination of graph-based and ridge penalty terms
#' 
#' @description 
#' Graph-constrained regression with penalty term being a linear combination of graph-based and ridge penalty terms.
#' 
#' See *Details* for model description and optimization problem formulation.  
#' 
#' @details 
#' Estimates coefficients of linear model of the formula: 
#' \deqn{y =  X\beta + Zb + \varepsilon}
#' where: 
#' - \eqn{y} - response,
#' - \eqn{X} - data matrix,
#' - \eqn{Z} - data matrix,
#' - \eqn{\beta} - regression coefficients, *not penalized* in estimation process,
#' - \eqn{b} - regression coefficients, *penalized* in estimation process and for whom there is, possibly a prior graph of similarity / graph of connections available.
#' 
#' The method uses a penalty being a linear combination of a graph-based and ridge penalty terms: 
#' \deqn{\beta_{est}, b_{est}= arg \; min_{\beta,b} \{ (y - X\beta - Zb)^T(y - X\beta - Zb) + \lambda_Qb^TQb +  \lambda_Rb^Tb \}}
#' where: 
#' - \eqn{Q} - a graph-originated penalty matrix; typically: a graph Laplacian matrix,
#' - \eqn{\lambda_Q} - regularization parameter for a graph-based penalty term
#' - \eqn{\lambda_R} - regularization parameter for ridge penalty term
#' 
#' The two regularization parameters, \eqn{\lambda_Q} and \eqn{\lambda_R}, are estimated as ML estimators from equivalent
#' Linear Mixed Model optimizaton problem formulation (see: References). 
#' 
#' * Graph-originated penalty term allows imposing similarity between coefficients based on graph information given.
#' * Ridge-originated penalty term facilitates parameters estimation: it reduces computational issues 
#'   arising from singularity in a graph-originated
#'   penalty matrix and yields plausible results in situations when graph information
#'   is not informative. 
#' 
#' Bootstrap confidence intervals computation is available (not set as a default option). 
#' 
#' @param Q graph-originated penalty matrix \eqn{(p \times p)}; typically: a graph Laplacian matrix
#' @param y response values matrix \eqn{(n \times 1)}
#' @param Z design matrix \eqn{(n \times p)} modeled as random effects variables (to be penalized in regression modeling); 
#' **assumed to be already standarized** 
#' @param X design matrix \eqn{(n \times k)} modeled as fixed effects variables (not to be penalized in regression modeling); 
#' if does not contain columns of 1s, such column will be added to be treated as intercept in a model 
#' 
#' @param optim.metod optimization method used to optimize \eqn{\lambda = (\lambda_Q, \lambda_R)}
#' * "rootSolve" (default) - optimizes by finding roots of non-linear equations by the Newton-Raphson method; from \code{rootSolve} package
#' * "sbplx" -  optimizes with the use of Subplex Algorithm: 'Subplex is a variant of Nelder-Mead that uses Nelder-Mead on a sequence of subspaces'; from \code{nloptr} package
#' @param rootSolve.x0 vector containing initial guesses for \eqn{\lambda = (\lambda_Q, \lambda_R)} used in "rootSolve" algorithm
#' @param rootSolve.Q0.x0 vector containing initial guess for \eqn{\lambda_R} used in "rootSolve" algorithm
#' @param sbplx.x0 vector containing initial guesses for \eqn{\lambda = (\lambda_Q, \lambda_R)} used in "sbplx" algorithm
#' @param sbplx.lambda.lo vector containing minimum values of \eqn{\lambda = (\lambda_Q, \lambda_R)} grid search in "sbplx" algorithm
#' @param sbplx.lambda.up vector containing mximum values of \eqn{\lambda = (\lambda_Q, \lambda_R)} grid search in "sbplx" algorithm
#' 
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
#' \item{lambda.R}{\eqn{\lambda_R} regularization parameter value}
#' \item{lambda.2}{\code{lambda.R}/\code{lambda.Q} value}
#' \item{boot.CI}{data frame with two columns, \code{lower} and \code{upper}, containing, respectively, values of lower and upper bootstrap confidence intervals for \eqn{b} regression coefficient estimates}
#' \item{obj.fn.val}{optimization problem objective function value}
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
#' riPEER.out <- riPEER(Q, y, Z, X)
#' plt.df <- data.frame(x = 1:p, y = riPEER.out$b.est)
#' ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line() + labs("b estimates")
#' }
#' 
#' \dontrun{
#' # riPEER with 0.95 bootstrap confidence intervals computation
#' riPEER.out <- riPEER(Q, y, Z, X, compute.boot.CI = TRUE, boot.R = 500)
#' plt.df <- data.frame(x = 1:p, 
#'                      y = riPEER.out$b.est, 
#'                      lo = riPEER.out$boot.CI[,1], 
#'                      up =  riPEER.out$boot.CI[,2])
#' ggplot(plt.df, aes(x = x, y = y, group = 1)) + geom_line() +  
#'   geom_ribbon(aes(ymin=lo, ymax=up), alpha = 0.3)
#' }
#' 
#' @references 
#' Karas, M., Brzyski, D., Dzemidzic, M., J., Kareken, D.A., Randolph, T.W., Harezlak, J. (2017). 
#' Brain connectivity-informed regularization methods for regression. doi: https://doi.org/10.1101/117945 
#' 
#' @import boot
#' @import nloptr
#' @import rootSolve
#' @export
#' 
riPEER <- function(Q, y, Z, X = NULL,
                   optim.metod = "rootSolve",
                   rootSolve.x0 = c(0.00001, 0.00001),
                   rootSolve.Q0.x0 = 0.00001,
                   sbplx.x0 = c(1,1),
                   sbplx.lambda.lo = c(10^(-5), 10^(-5)),
                   sbplx.lambda.up = c(1e6, 1e6),
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
  
  # Add intercept to X design matrix 
  if (is.null(X)){
    X <- matrix(rep(1, n), ncol = 1)
    warning("Adding colums of 1 to X data matrix. Treated as intercept in a model.")
  } else {
    X.var.variance   <- apply(X, 2, stats::var)
    X.var.variance.0 <- which(X.var.variance < 1e-10)
    if (length(X.var.variance.0) == 1){
      if (verbose) message(paste0("1 variable with variance=0 already in X data. ", 
                                  "Treated as intercept in a model."))
    } else if (length(X.var.variance.0) == 0){
      warning(paste0("No variable with variance=0 already in X data. ", 
                     "Adding colums of 1 to X matrix to be treated as intercept in a model."))
      X <- cbind(1, X)
    } else {
      stop(paste0(length(X.var.variance.0), " variables with variance=0 present in X data. ", 
                  "Cannot have more than 1 intercept-like column in X."))
    }
  }
  
  # Regress out fixed effect coefficients 
  X.hash <- diag(dim(X)[1]) - X %*% solve(t(X) %*% X) %*% t(X)
  y.wave <- X.hash %*% y
  Z.wave <- X.hash %*% Z
  
  # Transform Laplacian-constrained problem into generalized Ridge-constrained problem 
  Q.svd <- svd(Q)
  S <- diag(Q.svd$d)
  V <- Q.svd$u
  Z.waveV <- Z.wave %*% V

  # Optimize  
  if(verbose) message(paste0("Running optimizer: ", optim.metod, "..."))
  opt.out <- riPEER.opt.lambda(S = S, y = y.wave, Z = Z.waveV, 
                               optim.metod = optim.metod, 
                               sbplx.x0, sbplx.lambda.lo, sbplx.lambda.up,
                               rootSolve.x0, rootSolve.Q0.x0,
                               verbose)
  opt.obj.fn.val <- opt.out$obj.fn.val
  opt.lambda.Q   <- opt.out$lambda[1]
  opt.lambda.R   <- opt.out$lambda[2]
  if(verbose) message(paste0("Model lambdas optimized.",
                             "\nobj func val: ", round(opt.obj.fn.val, 3), 
                             "\nlambda_Q: ", round(opt.lambda.Q, 3), 
                             "\nlambda_R: ", round(opt.lambda.R, 3)))
  
  # Derive b, beta coefficients
  opt.D <- opt.lambda.Q * S + opt.lambda.R * diag(p)
  b.est.gR <- as.vector(solve(t(Z.waveV) %*% Z.waveV + 1 * opt.D) %*% t(Z.waveV) %*% y.wave)
  b.est <- V %*% b.est.gR
  beta.est <- solve(t(X) %*% X) %*% t(X) %*% (y - Z %*% b.est)
  
  # bootstrap confidence intervals
  if (compute.boot.CI){
    message(paste0("riPEER: computing bootstrap CI with ", boot.R, " replications..."))
    b.est.boot <- function(data, indices, opt.D, V) {
      data.boot <- data[indices, ]
      y.wave.boot   <- as.matrix(data.boot[,1])
      Z.waveV.boot  <- as.matrix(data.boot[,2:(ncol(data))])
      
      b.est.gR.boot <- as.vector(solve(t(Z.waveV.boot) %*% Z.waveV.boot + 
                                         1 * opt.D) %*% t(Z.waveV.boot) %*% y.wave.boot)
      b.est.boot <- V %*% b.est.gR.boot
      return(as.vector(b.est.boot))
    } 
    if (boot.set.seed) set.seed(1)
    boot.out <- boot(data = cbind(y.wave, Z.waveV), statistic = b.est.boot,  
                     R = boot.R, opt.D = opt.D, V = V, parallel = "multicore", ncpus = 7)
    boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx){
      boot.ci.out <- (boot.ci(boot.out, type = "bca", index = idx, conf = boot.conf))$bca
      return(boot.ci.out[c(4,5)])
    })
    boot.CI <- do.call(rbind.data.frame, boot.CI.l)
    names(boot.CI) <- c("lower", "upper")
  } else {
    boot.CI <- NULL
  }
  
  res.list <- list(b.est = as.vector(b.est), 
                   beta.est = as.vector(beta.est),
                   lambda.Q = opt.lambda.Q,
                   lambda.R = opt.lambda.R,
                   lambda.2 = opt.lambda.R / opt.lambda.Q,
                   boot.CI = boot.CI,
                   obj.fn.val = opt.obj.fn.val)
  return(res.list)
}