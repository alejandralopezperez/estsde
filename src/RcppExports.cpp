// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// GMM_CKLS
Rcpp::List GMM_CKLS(Rcpp::NumericVector r, double delta, Rcpp::NumericVector guess, int maxiter, double tol1, double tol2);
RcppExport SEXP _estsde_GMM_CKLS(SEXP rSEXP, SEXP deltaSEXP, SEXP guessSEXP, SEXP maxiterSEXP, SEXP tol1SEXP, SEXP tol2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type guess(guessSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol1(tol1SEXP);
    Rcpp::traits::input_parameter< double >::type tol2(tol2SEXP);
    rcpp_result_gen = Rcpp::wrap(GMM_CKLS(r, delta, guess, maxiter, tol1, tol2));
    return rcpp_result_gen;
END_RCPP
}
// MCMC_CKLS
Rcpp::List MCMC_CKLS(const Rcpp::NumericVector& obs, int n_obs, int m_data, int niter, int total_iter, double alpha, double kappa, double sigma, double gamma);
RcppExport SEXP _estsde_MCMC_CKLS(SEXP obsSEXP, SEXP n_obsSEXP, SEXP m_dataSEXP, SEXP niterSEXP, SEXP total_iterSEXP, SEXP alphaSEXP, SEXP kappaSEXP, SEXP sigmaSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< int >::type m_data(m_dataSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type total_iter(total_iterSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC_CKLS(obs, n_obs, m_data, niter, total_iter, alpha, kappa, sigma, gamma));
    return rcpp_result_gen;
END_RCPP
}
// MCMC_ou
Rcpp::List MCMC_ou(const Rcpp::NumericVector& obs, int n_obs, int m_data, int niter, int total_iter, double alpha, double kappa, double sigma);
RcppExport SEXP _estsde_MCMC_ou(SEXP obsSEXP, SEXP n_obsSEXP, SEXP m_dataSEXP, SEXP niterSEXP, SEXP total_iterSEXP, SEXP alphaSEXP, SEXP kappaSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< int >::type m_data(m_dataSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type total_iter(total_iterSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(MCMC_ou(obs, n_obs, m_data, niter, total_iter, alpha, kappa, sigma));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP FiltroKalmanCKLS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP FiltroKalmanVAS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP simCKLS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP simVAS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_estsde_GMM_CKLS", (DL_FUNC) &_estsde_GMM_CKLS, 6},
    {"_estsde_MCMC_CKLS", (DL_FUNC) &_estsde_MCMC_CKLS, 9},
    {"_estsde_MCMC_ou", (DL_FUNC) &_estsde_MCMC_ou, 8},
    {"FiltroKalmanCKLS", (DL_FUNC) &FiltroKalmanCKLS, 17},
    {"FiltroKalmanVAS",  (DL_FUNC) &FiltroKalmanVAS,  16},
    {"simCKLS",          (DL_FUNC) &simCKLS,           7},
    {"simVAS",           (DL_FUNC) &simVAS,            6},
    {NULL, NULL, 0}
};

RcppExport void R_init_estsde(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
