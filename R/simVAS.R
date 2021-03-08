
#' Simulation of the Vasicek model
#'
#' \code{rVAS} simulates a single path of the Vasicek process using its known conditional
#' distribution. The parametric form of the Vasicek model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma dW_t.}
#'
#' @param n an integer indicating the sample size.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param X0 a single numeric, initial value of the process.
#' @param alpha a single numeric, drift parameter.
#' @param kappa a single numeric, rate of mean reversion.
#' @param sigma a single numeric, volatility parameter.
#'
#' @return A \code{ts} object, a numeric vector with the sample path of the SDE.
#'
#' @examples
#' x <- rVAS(360, Delta = 1/12, X0 = 0.09, alpha = 0.08, kappa = 0.9, sigma = 0.1)
#' plot(x)
#'
#' @export
#' @useDynLib estsde simVAS
rVAS <- function(n, Delta, X0, alpha, kappa, sigma){
  variance <- sigma^2 * ((1 - exp(-2 * kappa * Delta)) / (2 * kappa))
  r <- .Call("simVAS", as.double(X0), as.integer(n), as.double(alpha),
             as.double(kappa), as.double(variance), as.double(Delta))
  r <- ts(r, start = 0, deltat = Delta, names = "X")
  return(r)
}
