#include <R.h>
#include <Rinternals.h>
SEXP FiltroKalmanVAS(SEXP y, SEXP rti, SEXP delta, SEXP alpha, SEXP beta, SEXP A, SEXP Phi, SEXP R, 
                      SEXP Q, SEXP xp, SEXP Pp, SEXP xf, SEXP Pf, SEXP K, SEXP innov, SEXP sig){
  int i, n;
  double siginv, like, *py, *prti, *pxp, *pPp, *pxf, *pPf, *pK, *pinnov, *psig;
  double rdelta, ralpha, rbeta, rA, rPhi, rR, rQ;
  SEXP out, list, list_names;   
  char *names[8] = {"xp", "Pp", "xf", "Pf", "like", "innov", "sig", "Kn"};
  n = length(y);

  rdelta = *REAL(delta);
  ralpha = *REAL(alpha);
  rbeta  = *REAL(beta);
  rA     = *REAL(A);
  rPhi   = *REAL(Phi);
  rR     = *REAL(R);
  rQ     = *REAL(Q);

  PROTECT(out = allocVector(REALSXP, 1));
  py     = REAL(y);
  prti   = REAL(rti);
  pxp    = REAL(xp);
  pPp    = REAL(Pp);
  pxf    = REAL(xf);
  pPf    = REAL(Pf);
  pK     = REAL(K);
  pinnov = REAL(innov);
  psig   = REAL(sig);

  siginv = 1 / psig[0];
  like = log(psig[0]) + pinnov[0] * siginv * pinnov[0];
  for (i = 1; i < n; i++) {
    pxp[i]    = rPhi * pxf[i - 1] + ralpha - rbeta * prti[i];
    pPp[i]    = rPhi * pPf[i - 1] * rPhi + rQ / rdelta;
    psig[i]   = rA * pPp[i] * rA + rR;
    siginv    = 1 / psig[i];
    pK[i]     = pPp[i] * rA * siginv;
    pinnov[i] = py[i] - rA * pxp[i];
    pxf[i]    = pxp[i] + pK[i] * pinnov[i];
    pPf[i]    = pPp[i] - pK[i] * rA * pPp[i];
    like += log(psig[i]) + pinnov[i] * siginv * pinnov[i];
  }
  like = 0.5 * like;
  REAL(out)[0] = like;
  // nombres de la lista
  PROTECT(list_names = allocVector(STRSXP,8));
  for(i = 0; i < 8; i++)   
     SET_STRING_ELT(list_names,i,mkChar(names[i]));
  // añadimos los objetos a la lista
  PROTECT(list = allocVector(VECSXP, 8)); 
  SET_VECTOR_ELT(list, 0, xp);
  SET_VECTOR_ELT(list, 1, Pp);
  SET_VECTOR_ELT(list, 2, xf);
  SET_VECTOR_ELT(list, 3, Pf);
  SET_VECTOR_ELT(list, 4, out);
  SET_VECTOR_ELT(list, 5, innov);
  SET_VECTOR_ELT(list, 6, sig);
  SET_VECTOR_ELT(list, 7, K);
  // atributos de la lista
  setAttrib(list, R_NamesSymbol, list_names);
  UNPROTECT(3);
  
  return list;
}
