
#' MCMC estimation for the Vasicek model
#'
#' Parametric estimation for the Vasicek model using Markov Chain Monte Carlo and
#' involving data augmentation, as proposed in Elerian et al. (2001) and Eraker (2001).
#' The parametric form of the Vasicek model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension three indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model as an initial guess.
#' @param niter an integer, number of iterations.
#' @param burn_in an integer indicating the number of initial iterations to be discarded.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(987)
#' x <- rVAS(480, 1/12, 0, 0.08, 0.9, 0.1)
#' est.VAS.MCMC(x)
#' }
#'
#' @references
#' Elerian, O., Chib, S., and Shephard, N. (2001). Likelihood inference for discretely
#' observed nonlinear diffusions. Econometrica, 69(4):959–993.
#'
#' Eraker, B. (2001). MCMC analysis of diffusion models with application to finance.
#' Journal of Business & Economic Statistics, 19(2):177–191.
est.VAS.MCMC <- function(X, Delta = deltat(X), par = NULL, niter = 2000, burn_in = 500) {

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  est       <- MCMC_ou(X, n_obs = length(X) - 1, m_data = 5, niter = niter, total_iter = niter + burn_in,
                       alpha = par[1]*Delta, kappa = par[2]*Delta, sigma = par[3]*sqrt(Delta))
  theta.hat <- est$param / c(Delta, Delta, sqrt(Delta))
  se        <- est$error/ c(Delta, Delta, sqrt(Delta))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}


#' MCMC estimation for the CKLS model
#'
#' Parametric estimation for the CKLS model using Markov Chain Monte Carlo and
#' involving data augmentation, as proposed in Elerian et al. (2001) and Eraker (2001).
#' The parametric form of the CKLS model used here is given by
#' \deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma X_t^\gamma dW_t.}
#'
#' @param X a numeric vector, the sample path of the SDE.
#' @param Delta a single numeric, the time step between two consecutive observations.
#' @param par a numeric vector with dimension four indicating initial values of the
#' parameters. Defaults to NULL, fits a linear model using generalized least squares
#' with AR1 correlation and a power variance heteroscedasticity structure.
#' @param niter an integer, number of iterations.
#' @param burn_in an integer indicating the number of initial iterations to be discarded.
#'
#' @return A list containing a matrix with the estimated coefficients and the
#' associated standard errors.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(987)
#' x <- rCKLS(480, 1/12, 0.09, 0.08, 0.9, 1.2, 1.5)
#' est.CKLS.MCMC(x)
#' }
#'
#' @references
#' Elerian, O., Chib, S., and Shephard, N. (2001). Likelihood inference for discretely
#' observed nonlinear diffusions. Econometrica, 69(4):959–993.
#'
#' Eraker, B. (2001). MCMC analysis of diffusion models with application to finance.
#' Journal of Business & Economic Statistics, 19(2):177–191.
est.CKLS.MCMC <- function(X, Delta = deltat(X), par = NULL, niter = 4000, burn_in = 1000) {

  if (is.null(par)) {
    yt   <- diff(X) / Delta
    Xt   <- X[-length(X)]
    init <- nlme::gls(yt ~ Xt, correlation = nlme::corAR1(), weights = nlme::varPower(form = ~ Xt))
    par  <- c(init$coefficients[1], -init$coefficients[2], init$sigma*sqrt(Delta), init$modelStruct$varStruct[1])
  }

  est       <- MCMC_CKLS(X, n_obs = length(X) - 1, m_data = 5, niter = niter, total_iter = niter + burn_in,
                         alpha = par[1]*Delta, kappa = par[2]*Delta, sigma = par[3]*sqrt(Delta), gamma = par[4])
  theta.hat <- est$param[-5] / c(Delta, Delta, sqrt(Delta), 1)
  se        <- est$error[-5] / c(Delta, Delta, sqrt(Delta), 1)

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma", "gamma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estCEV"
  return(res)
}


