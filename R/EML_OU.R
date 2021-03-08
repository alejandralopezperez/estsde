
# VAS conditional density
dVAS <- function(X, Delta, X0, Theta, log = FALSE){
  media    <- Theta[1] / Theta[2] + (X0 - Theta[1] / Theta[2]) * exp(-Theta[2] * Delta)
  variance <- Theta[3]^2 * ((1 - exp(-2 * Theta[2] * Delta)) / (2 * Theta[2]))
  dnorm(X, mean = media, sd = sqrt(variance), log = log)
}

# - log like
VAS.lik <- function(X, Delta, Theta) {
  n <- length(X)
  -sum(dVAS(X = X[2:n], Delta = Delta, X0 = X[1:(n - 1)], Theta = Theta, log = TRUE))
}


#' ML estimation for the Vasicek model
#'
#' Parametric estimation for the Vasicek model using (exact) maximum likelihood.
#' The parametric form of the Vasicek model used here is given by
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
#' est.VAS.EML(x)
est.VAS.EML <- function(X, Delta = deltat(X), par = NULL) {

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  est       <- optim(par, VAS.lik, X = X, Delta = Delta, hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}



