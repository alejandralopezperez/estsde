#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
SEXP simCKLS(SEXP r0, SEXP n, SEXP alpha, SEXP beta, SEXP sigma, SEXP gamma, SEXP delta){
  double pr0, palpha, pbeta, psigma, pgamma, pdelta;
  double *pr;
  double drift, diffusion, diffderiv, error;
  int pn, i;
  SEXP r;

  pn     = *INTEGER(n);
  pr0    = *REAL(r0);
  palpha = *REAL(alpha);
  pbeta  = *REAL(beta);
  psigma = *REAL(sigma);
  pgamma = *REAL(gamma);
  pdelta = *REAL(delta);

  PROTECT(r = allocVector(REALSXP, pn));
  pr = REAL(r);
  
  pr[0] = pr0;
  GetRNGstate();
  for (i = 1; i < pn; ++i)
  {
    drift     = palpha - pbeta * pr[i-1];
    diffusion = psigma * pow(pr[i-1], pgamma);
    diffderiv = psigma * pgamma * pow(pr[i-1], pgamma - 1);
    error     = rnorm(0, sqrt(pdelta));
    pr[i]     = pr[i-1] + drift * pdelta + diffusion * error + 0.5 * diffusion * diffderiv * (error * error - pdelta); 
  }
  PutRNGstate();

  UNPROTECT(1);

  return r;
}

