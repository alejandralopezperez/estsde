
denHP <- function(Delta, x, ssd, mu0, mu1, mu2, mu3, mu4, mu5, mu6){
  mu02 <- mu0^2
  mu03 <- mu0^3
  mu04 <- mu0^4
  mu05 <- mu0^5
  mu06 <- mu0^6
  mu12 <- mu1^2
  mu13 <- mu1^3
  mu22 <- mu2^2

  # Polinomios de Hermite
  H0 <- function(x) return(1.0)
  H1 <- function(x) return(-x)
  H2 <- function(x) return(x^2 - 1)
  H3 <- function(x) return(-x^3 + 3*x)
  H4 <- function(x) return(x^4 -6*x^2 + 3)
  H5 <- function(x) return(-x^5 + 10*x^3 - 15*x)
  H6 <- function(x) return(x^6 -15*x^4 +45*x^2 -15)

  eta1 <- -mu0*sqrt(Delta) - (2*mu0*mu1+mu2)*(Delta^1.5)/4 -
    (4*mu0*mu12 + 4*mu02*mu2 + 6*mu1*mu2 + 4*mu0*mu3 + mu4)*(Delta^2.5)/24

  eta2 <- (mu02+mu1)*Delta/2 + (6*mu02*mu1 + 4*mu12 + 7*mu0*mu2 + 2*mu3) *
    (Delta^2)/12 +(28*mu02*mu12 + 28*mu02*mu3 + 16*mu13 + 16*mu03*mu2 +
                     88*mu0*mu1*mu2 + 21*mu22 + 32*mu1*mu3 + 16*mu0*mu4 + 3*mu5)*(Delta^3)/96

  eta3 <- -(mu03 + 3*mu0*mu1 + mu2)*(Delta^1.5)/6 -
    (12*mu03 * mu1 + 28*mu0*mu12 + 22*mu02*mu2 + 24*mu1*mu2 + 14*mu0*mu3 + 3*mu4)*
    (Delta^2.5)/48

  eta4 <- (mu04 + 6*mu02 * mu1 + 3*mu12 + 4*mu0*mu2 + mu3)*(Delta^2)/24 +
    (20*mu04*mu1 + 50*mu03*mu2 + 100*mu02*mu12 + 50*mu02*mu3 + 23*mu0*mu4 +
       180*mu0*mu1*mu2 + 40*mu13 + 34*mu22 + 52*mu1*mu3 + 4*mu5)*(Delta^3.0)/240

  eta5 <- -(mu05 + 10*mu03*mu1 + 15*mu0*mu12 + 10*mu02*mu2 + 10*mu1*mu2 + 5*mu0*mu3 + mu4)*
    (Delta^2.5)/120

  eta6 <- (mu06 + 15*mu04*mu1 + 15*mu13 + 20*mu03*mu2 + 15*mu0*mu3 +
           45*mu02*mu12 + 10*mu22 + 15*mu02*mu3 + 60*mu0*mu1*mu2
           + 6*mu0*mu4 + mu5)*(Delta^3)/720

  val <- dnorm(x, 0, 1, FALSE) *
    (1 + eta1*H1(x) + eta2*H2(x) + eta3*H3(x) + eta4*H4(x) + eta5*H5(x) + eta6*H6(x))/ssd

  return(val)
}

HPlik <- function(X, Delta, Theta, mU0, mU1, mU2, mU3, mU4, mU5, mU6, lamperti, S){
  n   <- length(X)
  sd  <- sqrt(Delta)
  u0  <- lamperti(X[-n], Theta)
  f   <- lamperti(X[-1], Theta)
  ssd <- S(X[-1], Theta)*sd
  mu0 <- mU0(u0, Theta)
  mu1 <- mU1(u0, Theta)
  mu2 <- mU2(u0, Theta)
  mu3 <- mU3(u0, Theta)
  mu4 <- mU4(u0, Theta)
  mu5 <- mU5(u0, Theta)
  mu6 <- mU6(u0, Theta)
  val <- suppressWarnings(
    -sum(log(denHP(Delta, (f-u0)/sd, ssd, mu0, mu1, mu2, mu3, mu4, mu5, mu6)), na.rm = TRUE))
  return(val)
}


#' ML estimation for the Vasicek model (Hermite polynomial expansion)
#'
#' Parametric estimation for the Vasicek model using maximum likelihood and the
#' discretized version of the model, obtained with the Ait-Sahalia Hermite polynomial
#' expansion method. The parametric form of the Vasicek model used here is given by
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
#' est.VAS.HP(x)
#'
#' @references
#' Ait-Sahalia, Y. (2002). Maximum likelihood estimation of discretely sampled
#' diffusions: A closed-form approximation approach. Econometrica, 70(1):223–262.
est.VAS.HP <- function(X, Delta = deltat(X), par = NULL) {

  mU0  <- function(x, theta) theta[1] / theta[3] - theta[2] * x
  mU1  <- function(x, theta) - theta[2]
  mU2  <- function(x, theta) 0
  mU3  <- function(x, theta) 0
  mU4  <- function(x, theta) 0
  mU5  <- function(x, theta) 0
  mU6  <- function(x, theta) 0
  S    <- function(x, theta) theta[3]
  lamperti <- function(x, theta) x / theta[3]

  if (is.null(par)) {
    y    <- diff(X) / Delta
    init <- summary(lm(y ~ X[-length(X)]))
    par  <- c(init$coefficients[1,1], -init$coefficients[2,1], init$sigma * sqrt(Delta))
  }

  est       <- optim(par, HPlik, X = X, Delta = Delta, mU0 = mU0, mU1 = mU1, mU2 = mU2, mU3 = mU3,
                     mU4 = mU4, mU5 = mU5, mU6 = mU6, lamperti = lamperti, S = S, hessian = TRUE)
  theta.hat <- est$par
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estVAS"
  return(res)
}


#' ML estimation for the CKLS model (Hermite polynomial expansion)
#'
#' Parametric estimation for the CKLS model using maximum likelihood and the
#' discretized version of the model, obtained with the Ait-Sahalia Hermite polynomial
#' expansion method. The parametric form of the CKLS model used here is given by
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
#' est.CKLS.HP(x)
#'
#' @references
#' Ait-Sahalia, Y. (2002). Maximum likelihood estimation of discretely sampled
#' diffusions: A closed-form approximation approach. Econometrica, 70(1):223–262.
est.CKLS.HP <- function(X, Delta = deltat(X), par = NULL) {

  mU0 <- function(x, theta) {
    theta[4]/(2*(theta[4] - 1)*x) - theta[2]*(theta[4]-1)*x + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * x^(theta[4]/(theta[4]-1))
  }
  mU1 <- function(x, theta) {
    -theta[4]/(2*(theta[4] - 1)*x^2) - theta[2]*(theta[4]-1) +
      theta[1]*theta[3]^(1/(theta[4]-1)) * (theta[4]-1)^(theta[4]/(theta[4]-1)) *
      (theta[4]/(theta[4]-1)) * x^(theta[4]/(theta[4]-1) - 1)
  }
  mU2 <- function(x, theta) {
    theta[4]/((theta[4] - 1)*x^3) + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * (theta[4]/(theta[4]-1)^2) * x^(theta[4]/(theta[4]-1) - 2)
  }
  mU3 <- function(x, theta) {
    -3*theta[4]/((theta[4] - 1)*x^4) + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * (theta[4]*(2-theta[4])/(theta[4]-1)^3) *
      x^(theta[4]/(theta[4]-1) - 3)
  }
  mU4 <- function(x, theta) {
    12*theta[4]/((theta[4] - 1)*x^5) + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * (theta[4]*(2*theta[4]^2 - 7*theta[4] + 6)/(theta[4]-1)^4) *
      x^(theta[4]/(theta[4]-1) - 4)
  }
  mU5 <- function(x, theta) {
    -60*theta[4]/((theta[4] - 1)*x^6) + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * (theta[4]*(-6*theta[4]^3 + 29*theta[4]^2 -
       46*theta[4] + 24)/(theta[4]-1)^5) * x^(theta[4]/(theta[4]-1) - 5)
  }
  mU6 <- function(x, theta) {
    360*theta[4]/((theta[4] - 1)*x^7) + theta[1]*theta[3]^(1/(theta[4]-1)) *
      (theta[4]-1)^(theta[4]/(theta[4]-1)) * (theta[4]*(24*theta[4]^4 - 146*theta[4]^3 +
       329*theta[4]^2 - 326*theta[4] + 120)/(theta[4]-1)^6) * x^(theta[4]/(theta[4]-1) - 6)
  }
  S   <- function(x, theta) theta[3] * x^theta[4]
  lamperti <- function(x, theta) x^(1 - theta[4]) / (theta[3] * (theta[4] - 1))


  if (is.null(par)) {
    yt   <- diff(X) / Delta
    Xt   <- X[-length(X)]
    init <- nlme::gls(yt ~ Xt, correlation = nlme::corAR1(), weights = nlme::varPower(form = ~ Xt))
    par  <- c(init$coefficients[1], -init$coefficients[2], init$sigma*sqrt(Delta), init$modelStruct$varStruct[1])
  }

  est       <- optim(par, HPlik, X = X, Delta = Delta, mU0 = mU0, mU1 = mU1, mU2 = mU2, mU3 = mU3,
                     mU4 = mU4, mU5 = mU5, mU6 = mU6, lamperti = lamperti, S = S, hessian = TRUE)
  theta.hat <- abs(est$par)
  se        <- sqrt(diag(solve(est$hessian)))

  coeff <- cbind(theta.hat, se)
  rownames(coeff) <- c("alpha", "kappa", "sigma", "gamma")
  colnames(coeff) <- c("Estimate", "Std. Error")
  res <- list(coefficients = coeff)
  class(res) <- "estCEV"
  return(res)
}
