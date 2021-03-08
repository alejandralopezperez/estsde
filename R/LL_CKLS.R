
# density LL -------
denLL_CKLS <- function(X, X0, Delta, Theta, d, dx, log = FALSE) {
  drif.deriv <- dx(X0, Theta)
  drift      <- d(X0, Theta)
  K          <- log(1 + drift * (exp(drif.deriv * Delta) - 1)/(X0 * drif.deriv))/Delta
  esperanza  <- X0 + drift/drif.deriv * (exp(drif.deriv * Delta) - 1)
  variance   <- (exp(2 * K * Delta) - 1)/(2 * K)
  dnorm(X, mean = esperanza, sd = sqrt(variance), log = log)
}

# likelihood LL ----
LLlik_CKLS <- function(X, Delta, Theta, d, dx) {
  U <- X^(1 - Theta[4]) / (Theta[3] * (Theta[4] - 1))  # Lamperti
  n <- length(X)
  suppressWarnings(
    -sum(denLL_CKLS(X = U[2:n], X0 = U[1:(n-1)], Delta = Delta, Theta = Theta,
                    d = d, dx = dx, log = TRUE)  + log(1/(Theta[3] * X[2:n]^Theta[4])),
         na.rm = TRUE))
}


#' ML estimation for the CKLS model (local linearization)
#'
#' Parametric estimation for the CKLS model using maximum likelihood and
#' the discretized version of the model, obtained with the local linearization
#' method. The parametric form of the CKLS model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma X_t^\gamma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension four indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model using generalized least squares
#' with AR1 correlation and a power variance heteroscedasticity structure.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' x <- rCKLS(360, 1/12, 0.09, 0.08, 0.9, 1.2, 1.5)
#' est.CKLS.LL(x)
#'
#' @references
#' Ozaki, T. (1992). A bridge between nonlinear time series models and
#' nonlinear stochastic dynamical systems: a local linearization approach.
#' Statistica Sinica, pages 113–135.
#'
#' Shoji, I. and Ozaki, T. (1998). Estimation for nonlinear stochastic
#' differential equations by a local linearization method. Stochastic
#' Analysis and Applications, 16(4):733–752.
est.CKLS.LL <- function(X, Delta = deltat(X), par = NULL) {

  if (is.null(par)) {
    yt   <- diff(X) / Delta
    Xt   <- X[-length(X)]
    init <- nlme::gls(yt ~ Xt, correlation = nlme::corAR1(), weights = nlme::varPower(form = ~ Xt))
    par  <- c(init$coefficients[1], -init$coefficients[2], init$sigma*sqrt(Delta), init$modelStruct$varStruct[1])
  }

  d  <- function(x, theta) {
    theta[4]/(2*(theta[4] - 1)*x) - theta[2]*(theta[4]-1)*x + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * x^(theta[4]/(theta[4]-1))
  }
  dx <- function(x, theta) {
    -theta[4]/(2*(theta[4] - 1)*x^2) - theta[2]*(theta[4]-1) +
      theta[1]*theta[3]^(1/(theta[4]-1)) * (theta[4]-1)^(theta[4]/(theta[4]-1)) *
      (theta[4]/(theta[4]-1)) * x^(theta[4]/(theta[4]-1) - 1)
  }

  est       <- optim(par, LLlik_CKLS, X = X, Delta = Delta, d = d, dx = dx, method = 'BFGS', hessian = TRUE)
  theta.hat <- abs(est$par)
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma", "gamma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estCEV"
  return(res)
}

