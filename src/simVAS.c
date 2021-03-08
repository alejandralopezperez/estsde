#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
SEXP simVAS(SEXP r0, SEXP n, SEXP alpha, SEXP beta, SEXP varianza, SEXP delta){
  double pr0, palpha, pbeta, pvarianza, pdelta;
  double *pr;
  double media;
  int pn, i;
  SEXP r;

  pn        = *INTEGER(n);
  pr0       = *REAL(r0);
  palpha    = *REAL(alpha);
  pbeta     = *REAL(beta);
  pvarianza = *REAL(varianza);
  pdelta    = *REAL(delta);

  PROTECT(r = allocVector(REALSXP, pn));
  pr = REAL(r);
  
  pr[0] = pr0;
  GetRNGstate();
  for (i = 1; i < pn; ++i)
  {
    media = palpha / pbeta + (pr[i-1] - palpha / pbeta) * exp(-pbeta * pdelta);
    pr[i] = rnorm(media, sqrt(pvarianza));
  }
  PutRNGstate();

  UNPROTECT(1);

  return r;
}
