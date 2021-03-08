
#' @useDynLib estsde FiltroKalmanCKLS
FK.CKLS <- function (y, rti, delta, alpha, beta, gamma, mu0, Sigma0, cQ) {
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
  Pp[1]    <- Phi * P00 * t(Phi) + Q * (rti[1]^(2 * gamma)) / delta
  sig[1]   <- A * Pp[1] * t(A) + R
  siginv   <- 1 / sig[1]
  K[1]     <- Pp[1] * t(A) * siginv
  innov[1] <- y[1] - A * xp[1]
  xf[1]    <- xp[1] + K[1] * innov[1]
  Pf[1]    <- Pp[1] - K[1] * A * Pp[1]
  like     <- log(sig[1]) + t(innov[1]) * siginv * innov[1]

  out <- .Call("FiltroKalmanCKLS", as.double(y), as.double(rti), as.double(delta), as.double(alpha), as.double(beta),
               as.double(gamma), as.double(A), as.double(Phi), as.double(R), as.double(Q), as.double(xp),
               as.double(Pp), as.double(xf), as.double(Pf), as.double(K), as.double(innov), as.double(sig))
  return(out)
}


KFlik_CKLS <- function(X, Delta, Theta, mu0, Sigma0) {
  y     <- diff(X) / Delta
  num   <- length(y)
  alpha <- Theta[1]
  beta  <- Theta[2]
  cQ    <- Theta[3]
  gamm  <- Theta[4]

  est <- FK.CKLS(y = y, rti = X, delta = Delta, alpha = alpha, beta = beta,
                 gamma = gamm, mu0 = mu0, Sigma0 = Sigma0, cQ = cQ)
  return(est$like)
}

#' ML estimation for the CKLS model (Kalman filter)
#'
#' Parametric estimation for the CKLS model using the Kalman filter algorithm, where
#' the discretized version of the model is obtained with the Euler-Maruyama method.
#' The parametric form of the CKLS model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma X_t^\gamma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension four indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model using generalized least squares
#' with AR1 correlation and a power variance heteroscedasticity structure.
#' @param mu0 a single numeric, the initial mean. Defaults to zero.
#' @param Sigma0 a single numeric, the initial variance. Defaults to one.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' x <- rCKLS(360, 1/12, 0.09, 0.08, 0.9, 1.2, 1.5)
#' est.CKLS.KF(x)
est.CKLS.KF <- function(X, Delta = deltat(X), par = NULL, mu0 = 0, Sigma0 = 1) {

  if (is.null(par)) {
    yt   <- diff(X) / Delta
    Xt   <- X[-length(X)]
    init <- nlme::gls(yt ~ Xt, correlation = nlme::corAR1(), weights = nlme::varPower(form = ~ Xt))
    par  <- c(init$coefficients[1], -init$coefficients[2], init$sigma*sqrt(Delta), init$modelStruct$varStruct[1])
  }

  est       <- optim(par, KFlik_CKLS, X = X, Delta = Delta, mu0 = mu0, Sigma0 = Sigma0,
                     method = 'BFGS', hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma", "gamma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estCEV"
  return(res)
}
