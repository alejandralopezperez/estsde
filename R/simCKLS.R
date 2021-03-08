
#' Simulation of the CKLS model
#'
#' \code{rCKLS} simulates a single path of the CKLS process using the Milstein scheme.
#' The parametric form of the CKLS model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma X_t^\gamma dW_t.}
#'
#' @param n an integer indicating the sample size.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param X0 a single numeric, initial value of the process.
#' @param alpha a single numeric, drift parameter.
#' @param kappa a single numeric, rate of mean reversion.
#' @param sigma a single numeric, volatility parameter.
#' @param gamma a single numeric, proportional volatility exponent.
#'
#' @return A \code{ts} object, a numeric vector with the sample path of the SDE.
#'
#' @examples
#' x <- rCKLS(360, Delta = 1/12, X0 = 0.09, alpha = 0.08, kappa = 0.9, sigma = 1.2, gamma = 1.5)
#' plot(x)
#'
#' @references
#' Chan, K. C., Karolyi, G. A., Longstaff, F. A., and Sanders, A. B. (1992).
#' An empirical comparison of alternative models of the short-term interest rate.
#' The journal of finance, 47(3):1209–1227.
#'
#' @export
#' @useDynLib estsde simCKLS
rCKLS <- function(n, Delta, X0, alpha, kappa, sigma, gamma){
  r <- .Call("simCKLS", as.double(X0), as.integer(n), as.double(alpha), as.double(kappa),
             as.double(sigma), as.double(gamma), as.double(Delta))
  r <- ts(r, start = 0, deltat = Delta, names = "X")
  return(r)
}
