# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

GMM_CKLS <- function(r, delta, guess, maxiter, tol1, tol2) {
    .Call('_estsde_GMM_CKLS', PACKAGE = 'estsde', r, delta, guess, maxiter, tol1, tol2)
}

MCMC_CKLS <- function(obs, n_obs, m_data, niter, total_iter, alpha, kappa, sigma, gamma) {
    .Call('_estsde_MCMC_CKLS', PACKAGE = 'estsde', obs, n_obs, m_data, niter, total_iter, alpha, kappa, sigma, gamma)
}

MCMC_ou <- function(obs, n_obs, m_data, niter, total_iter, alpha, kappa, sigma) {
    .Call('_estsde_MCMC_ou', PACKAGE = 'estsde', obs, n_obs, m_data, niter, total_iter, alpha, kappa, sigma)
}

