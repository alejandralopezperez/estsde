% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HP.R
\name{est.VAS.HP}
\alias{est.VAS.HP}
\title{ML estimation for the Vasicek model (Hermite polynomial expansion)}
\usage{
est.VAS.HP(X, Delta = deltat(X), par = NULL)
}
\arguments{
\item{X}{a numeric vector, the sample path of the SDE.}

\item{Delta}{a single numeric, the time step between two consecutive observations.}

\item{par}{a numeric vector with dimension three indicating initial values of the
parameters. Defaults to NULL, fits a linear model as an initial guess.}
}
\value{
A list containing a matrix with the estimated coefficients and the
associated standard errors.
}
\description{
Parametric estimation for the Vasicek model using maximum likelihood and the
discretized version of the model, obtained with the Ait-Sahalia Hermite polynomial
expansion method. The parametric form of the Vasicek model used here is given by
\deqn{dX_t = (\alpha - \kappa X_t)dt + \sigma dW_t.}
}
\examples{
x <- rVAS(360, 1/12, 0, 0.08, 0.9, 0.1)
est.VAS.HP(x)

}
\references{
Ait-Sahalia, Y. (2002). Maximum likelihood estimation of discretely sampled
diffusions: A closed-form approximation approach. Econometrica, 70(1):223–262.
}
