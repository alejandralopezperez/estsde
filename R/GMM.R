
#' Generalized method of moments estimator for the Vasicek model
#'
#' Parametric estimation for the Vasicek model using the Generalized Method
#' of Moments. The parametric form of the Vasicek model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension three indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model as an initial guess.
#' @param maxiter an integer, the maximum number of iterations.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' x <- rVAS(360, 1/12, 0, 0.08, 0.9, 0.1)
#' est.VAS.GMM(x)
#'
#' @references
#' Hansen, L. P. (1982). Large sample properties of generalized method of moments
#' estimators. Econometrica, pages 1029–1054.
est.VAS.GMM <- function(X, Delta = deltat(X), par = NULL, maxiter = 25) {

  tol1 <- 0.001
  tol2 <- 0.001

  ft <- function(x, y, Theta, Delta){
    c.mean <- Theta[1]/Theta[2] + (y-Theta[1]/Theta[2])*exp(-Theta[2]*Delta)
    c.var <- Theta[3]^2 * (1-exp(-2*Theta[2]*Delta))/(2*Theta[2])
    cbind(x-c.mean, y*(x-c.mean), c.var-(x-c.mean)^2, y*(c.var-(x-c.mean)^2))
  }

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  n  <- length(X)
  gn <- function(theta) apply(ft(X[2:n], X[1:(n - 1)], theta, Delta), 2, mean)
  Q  <- function(theta) sum(gn(theta)^2)
  S  <- function(j, theta) ((t(ft(X[(j + 2):n], X[(j + 1):(n - 1)], theta, Delta)) %*% ft(X[2:(n - j)], X[1:(n - j - 1)], theta, Delta))/n)
  Q1 <- function(theta, W) gn(theta) %*% W %*% gn(theta)

  ell <- n - 2
  w   <- 1 - (1:ell)/(ell + 1)

  theta0 <- optim(par, Q, method = "BFGS")$par

  go   <- TRUE
  iter <- 0
  while (go) {
    iter  <- iter + 1
    S.hat <- S(0, theta0)
    for (i in 1:ell) S.hat = S.hat + w[i] * (S(i, theta0) + t(S(i, theta0)))
    W      <- solve(S.hat)
    est    <- optim(theta0, Q1, W = W, method = "BFGS")
    theta1 <- est$par
    if (sum(abs(theta0 - theta1)) < tol1 || est$value < tol2 || iter > maxiter) go <- FALSE
    theta0 <- theta1
  }

  dhat <- numDeriv::jacobian(gn, theta0)
  se   <- sqrt(diag(solve(t(dhat) %*% W %*% dhat)) / n)

  coeff <- cbind(theta0, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}


#' Generalized method of moments estimator for the CKLS model
#'
#' Parametric estimation for the CKLS model using the Generalized Method
#' of Moments. The parametric form of the CKLS model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma X_t^\gamma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension four indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model using generalized least squares
#' with AR1 correlation and a power variance heteroscedasticity structure.
#' @param maxiter an integer, the maximum number of iterations.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' x <- rCKLS(360, 1/12, 0.09, 0.08, 0.9, 1.2, 1.5)
#' est.CKLS.GMM(x)
#'
#' @references
#' Hansen, L. P. (1982). Large sample properties of generalized method of moments
#' estimators. Econometrica, pages 1029–1054.
#'
#' Chan, K. C., Karolyi, G. A., Longstaff, F. A., and Sanders, A. B. (1992).
#' An empirical comparison of alternative models of the short-term interest rate.
#' The journal of finance, 47(3):1209–1227.
est.CKLS.GMM <- function(X, Delta = deltat(X), par = NULL, maxiter = 25) {

  tol1 <- 0.001
  tol2 <- 0.001

  ft <- function(x, y, Theta, Delta){
    c.mean <- y + (Theta[1] - Theta[2]*y) * Delta
    c.var <- (Theta[3]^2) * (y^(2 * Theta[4])) * Delta
    cbind(x-c.mean,y*(x-c.mean), c.var-(x-c.mean)^2, y*(c.var-(x-c.mean)^2))
  }

  if (is.null(par)) {
    yt   <- diff(X) / Delta
    Xt   <- X[-length(X)]
    init <- nlme::gls(yt ~ Xt, correlation = nlme::corAR1(), weights = nlme::varPower(form = ~ Xt))
    par  <- c(init$coefficients[1], -init$coefficients[2], init$sigma*sqrt(Delta), init$modelStruct$varStruct[1])
  }

  n  <- length(X)
  gn <- function(theta) apply(ft(X[2:n], X[1:(n - 1)], theta, Delta), 2, mean)

  est <- GMM_CKLS(X, Delta, par, maxiter, tol1, tol2)

  dhat <- numDeriv::jacobian(gn, est$par)
  se   <- sqrt(diag(solve(t(dhat) %*% est$W %*% dhat)) / n)

  coeff <- cbind(est$par, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma", "gamma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estCEV"
  return(res)
}

