
#' @useDynLib estsde FiltroKalmanVAS
KF.VAS <- function (y, rti, delta, alpha, beta, mu0, Sigma0, cQ) {
  A   <- 1
  Phi <- 0
  R   <- 0
  n   <- length(y)
  Q   <- t(cQ) %*% cQ

  xp    <- numeric(n)
  Pp    <- numeric(n)
  xf    <- numeric(n)
  Pf    <- numeric(n)
  K     <- numeric(n)
  innov <- numeric(n)
  sig   <- numeric(n)

  x00 <- mu0
  P00 <- Sigma0

  xp[1]    <- Phi * x00 + alpha - beta * rti[1]
  Pp[1]    <- Phi * P00 * t(Phi) + Q / delta
  sig[1]   <- A * Pp[1] * t(A) + R
  siginv   <- 1 / sig[1]
  K[1]     <- Pp[1] * t(A) * siginv
  innov[1] <- y[1] - A * xp[1]
  xf[1]    <- xp[1] + K[1] * innov[1]
  Pf[1]    <- Pp[1] - K[1] * A * Pp[1]
  like     <- log(sig[1]) + t(innov[1]) * siginv * innov[1]

  out <- .Call("FiltroKalmanVAS", as.double(y), as.double(rti), as.double(delta), as.double(alpha), as.double(beta),
               as.double(A), as.double(Phi), as.double(R), as.double(Q), as.double(xp), as.double(Pp),
               as.double(xf), as.double(Pf), as.double(K), as.double(innov), as.double(sig))
  return(out)
}


KFlik_OU <- function(X, Delta, Theta, mu0, Sigma0) {
  y     <- diff(X) / Delta
  num   <- length(y)
  alpha <- Theta[1]
  beta  <- Theta[2]
  cQ    <- Theta[3]

  est <- KF.VAS(y = y, rti = X, delta = Delta, alpha = alpha,
                beta = beta, mu0 = mu0, Sigma0 = Sigma0, cQ = cQ)
  return(est$like)
}


#' ML estimation for the Vasicek model (Kalman filter)
#'
#' Parametric estimation for the Vasicek model using the Kalman filter algorithm, where
#' the discretized version of the model is obtained with the Euler-Maruyama method.
#' The parametric form of the Vasicek model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension three indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model as an initial guess.
#' @param mu0 a single numeric, the initial mean. Defaults to zero.
#' @param Sigma0 a single numeric, the initial variance. Defaults to one.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' x <- rVAS(360, 1/12, 0, 0.08, 0.9, 0.1)
#' est.VAS.KF(x)
est.VAS.KF <- function(X, Delta = deltat(X), par = NULL, mu0 = 0, Sigma0 = 1) {

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  est       <- optim(par, KFlik_OU, X = X, Delta = Delta, mu0 = mu0, Sigma0 = Sigma0,
                     method = 'BFGS', hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}

