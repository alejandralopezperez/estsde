# estsde

This software provides different parametric estimation methods for
stochastic differential equations (SDE), namely the Vasicek and CKLS
models, as well as functions to simulate paths of the SDEs. The methods
included are:

  - Exact maximum likelihood (EML)
  - Euler method (DML)
  - Local linearization (LL)
  - Hermite polynomial expansion (HP)
  - Generalized Method of Moments (GMM)
  - Kalman Filter (KF)
  - Markov Chain Monte Carlo (MCMC)

## Installation

Get the package from GitHub:

``` r
# Install the package
library(devtools)
install_github("alejandralopezperez/estsde")
# Load package
library(estsde)
```

## Usage

The following is an example of data generation and estimation:

``` r
# Generate data from a CKLS model
set.seed(789)
x <- rCKLS(480, 1/12, 0.09, 0.08, 0.9, 1.2, 1.5)
# Parameter estimates
est.CKLS.KF(x)
#> 
#> Call: 
#> dX_t = (alpha - kappa X_t)dt + sigma X_t^gamma dW_t;
#> 
#> Coefficients: 
#>       Estimate Std. Error
#> alpha   0.0795     0.0192
#> kappa   0.8979     0.2530
#> sigma   1.3393     0.3890
#> gamma   1.5496     0.1176
```
