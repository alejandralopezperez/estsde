
# likelihood -------
DMLlik_OU <- function(X, Delta, Theta){
  suppressWarnings(
  -sum(dnorm(X[-1], X[-length(X)] + (Theta[1] - Theta[2]*X[-length(X)]) * Delta,
             sqrt(Delta) * Theta[3], log = TRUE)))
}

#' ML estimation for the Vasicek model (Euler method)
#'
#' Parametric estimation for the Vasicek model using maximum likelihood and
#' the discretized version of the model, obtained with the Euler-Maruyama
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
#' est.VAS.DML(x)
est.VAS.DML <- function(X, Delta = deltat(X), par = NULL) {

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  est       <- optim(par, DMLlik_OU, X = X, Delta = Delta, hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}

