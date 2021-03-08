
# density -------
denLL_OU <- function(X, X0, Delta, Theta, d, dx, log = FALSE) {
  drif.deriv <- dx(X0, Theta)
  drift      <- d(X0, Theta)
  K          <- log(1 + drift * (exp(drif.deriv * Delta) - 1)/(X0 * drif.deriv))/Delta
  esperanza  <- X0 + drift/drif.deriv * (exp(drif.deriv * Delta) - 1)
  variance   <- Theta[3]^2 * (exp(2 * K * Delta) - 1)/(2 * K)
  dnorm(X, mean = esperanza, sd = sqrt(variance), log = log)
}

# likelihood ----
LLlik_OU <- function(X, Delta, Theta, d, dx) {
  n <- length(X)
  suppressWarnings(
    -sum(denLL_OU(X = X[2:n], X0 = X[1:(n-1)], Delta = Delta, Theta = Theta,
                  d = d, dx = dx, log = TRUE), na.rm = TRUE))
}

#' ML estimation for the Vasicek model (local linearization)
#'
#' Parametric estimation for the Vasicek model using maximum likelihood and
#' the discretized version of the model, obtained with the local linearization
#' method. The parametric form of the Vasicek model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension three indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model as an initial guess.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' x <- rVAS(360, 1/12, 0, 0.08, 0.9, 0.1)
#' est.VAS.LL(x)
#'
#' @references
#' Ozaki, T. (1992). A bridge between nonlinear time series models and
#' nonlinear stochastic dynamical systems: a local linearization approach.
#' Statistica Sinica, pages 113–135.
#'
#' Shoji, I. and Ozaki, T. (1998). Estimation for nonlinear stochastic
#' differential equations by a local linearization method. Stochastic
#' Analysis and Applications, 16(4):733–752.
est.VAS.LL <- function(X, Delta = deltat(X), par = NULL) {

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  d  <- function(x, Theta) Theta[1] - Theta[2] * x
  dx <- function(x, Theta) - Theta[2]

  est       <- optim(par, LLlik_OU, X = X, Delta = Delta, d = d, dx = dx, hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}

