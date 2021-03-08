
# likelihood -------
DMLlik_CKLS <- function(X, Delta, Theta){
  -sum(dnorm(X[-1], X[-length(X)] + (Theta[1] - Theta[2]*X[-length(X)]) * Delta,
             sqrt(Delta) * Theta[3] * X[-length(X)]^Theta[4], log = TRUE))
}


#' ML estimation for the CKLS model (Euler method)
#'
#' Parametric estimation for the CKLS model using maximum likelihood and
#' the discretized version of the model, obtained with the Euler-Maruyama
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
#' est.CKLS.DML(x)
est.CKLS.DML <- function(X, Delta = deltat(X), par = NULL) {

  if (is.null(par)) {
    yt   <- diff(X) / Delta
    Xt   <- X[-length(X)]
    init <- nlme::gls(yt ~ Xt, correlation = nlme::corAR1(), weights = nlme::varPower(form = ~ Xt))
    par  <- c(init$coefficients[1], -init$coefficients[2], init$sigma*sqrt(Delta), init$modelStruct$varStruct[1])
  }

  est       <- optim(par, DMLlik_CKLS, X = X, Delta = Delta, hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma", "gamma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estCEV"
  return(res)
}


