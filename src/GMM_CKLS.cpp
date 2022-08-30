// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include <RcppArmadillo.h>
#include <nloptrAPI.h>

using namespace arma;
using namespace Rcpp;

typedef struct {
  Rcpp::List dato;
} my_data;

template <class T>
inline double do_mean( T& x ) {
  return mean(x);
}


Rcpp::NumericVector col_means( NumericMatrix X ) {

  int nCols = X.ncol();
  NumericVector out = no_init(nCols);

  for(int j=0; j < nCols; j++) {
    NumericMatrix::Column tmp = X(_, j);
    out[j] = do_mean( tmp );
  }

  return out;
}


Rcpp::NumericMatrix f(double theta1, double theta2, double theta3, double theta4,
                      double Delta, Rcpp::NumericVector rt, int j, int k) {

  int n = rt.size();
  Rcpp::NumericMatrix g(n-1-k-j, 4);
  Rcpp::NumericVector r(n-1-k-j), r1(n-1-k-j), ehat(n-1-k-j);

  r    = rt[Rcpp::Range(1+j, n-1-k)];
  r1   = rt[Rcpp::Range(0+j, n-2-k)];
  ehat = r - r1 - (theta1 - theta2 * r1) * Delta;

  g(_,0) = ehat;
  g(_,1) = ehat * r1;
  g(_,2) = pow(ehat, 2) - Delta * pow(theta3, 2) * pow(r1, 2*theta4);
  g(_,3) = g(_,2) * r1;

  return g;
}


Rcpp::NumericVector H(double theta1, double theta2, double theta3,
                      double theta4, double Delta, Rcpp::NumericVector rt) {

  return col_means(f(theta1, theta2, theta3, theta4, Delta, rt, 0, 0));
}


double Q(double theta1, double theta2, double theta3,
         double theta4, Rcpp::List delta_r) { // double Delta, Rcpp::NumericVector rt

  double Delta = delta_r[0];
  Rcpp::NumericVector rt = delta_r[1];
  return sum(pow(H(theta1, theta2, theta3, theta4, Delta, rt), 2));
}


double Q_w(unsigned n, const double *x, double *grad, void *my_func_data){

  my_data *temp = (my_data *) my_func_data;
  Rcpp::List dato = temp->dato;
  double result = Q(x[0], x[1], x[2], x[3], dato);
  return result;
}


arma::mat S(double theta1, double theta2, double theta3, double theta4,
            double Delta, Rcpp::NumericVector rt, int j, int k) {

  int n = rt.size();
  arma::mat tmp1 = Rcpp::as<arma::mat>(f(theta1, theta2, theta3, theta4, Delta, rt, j, 0));
  arma::mat tmp2 = Rcpp::as<arma::mat>(f(theta1, theta2, theta3, theta4, Delta, rt, 0, k));
  arma::mat res = (tmp1.t() * tmp2) / n; // 4x4
  return res;
}


double Q1(double theta1, double theta2, double theta3,
          double theta4, Rcpp::List delta_r_W) {  // double Delta, Rcpp::NumericVector rt, arma::mat W

  double Delta = delta_r_W[0];
  Rcpp::NumericVector rt = delta_r_W[1];
  arma::mat W = delta_r_W[2];
  arma::vec tmp = Rcpp::as<arma::vec>(H(theta1, theta2, theta3, theta4, Delta, rt));
  arma::mat res = tmp.t() * W * tmp; // matrix 1x1
  return res(0);
}


double Q1_w(unsigned n, const double *x, double *grad, void *my_func_data) {
  my_data *temp = (my_data *) my_func_data;
  Rcpp::List dato = temp->dato;
  double result = Q1(x[0], x[1], x[2], x[3], dato);
  return result;
}


Rcpp::List optQ(Rcpp::NumericVector x, double Delta, Rcpp::NumericVector r) {

  double lb[4] = {0, 0, 0, 0};                               /* lower bounds */
  double ub[4] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL};   /* upper bounds */
  nlopt_opt opt;
  opt = nlopt_create(NLOPT_LN_NELDERMEAD, 4);                /* algorithm and dimensionality */
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  Rcpp::List lista = Rcpp::List::create(Rcpp::Named("delta") = Delta,
                                        Rcpp::Named("r") = r);
  my_data &temp = (my_data &) lista;
  nlopt_set_min_objective(opt, Q_w, &temp);

  nlopt_set_xtol_rel(opt, 1e-8);
  double minf;                                               /* the minimum objective value, upon return */
  // nlopt_result result = nlopt_optimize(opt, x.begin(), &minf);
  nlopt_optimize(opt, x.begin(), &minf);
  nlopt_destroy(opt);                                        /* dispose of the nlopt_opt object */

  Rcpp::NumericVector theta(4);
  for (int i = 0; i < 4; ++i) theta[i] = x[i];

  double res = minf;
  // Rcpp::Rcout << "nlopt_result: " << result << "\n";
  // Rcpp::Rcout << "theta_hat: " << theta << "\n\n";

  return Rcpp::List::create(Rcpp::Named("par") = theta, Rcpp::Named("val") = res);
}

Rcpp::List optQ1(Rcpp::NumericVector x, double Delta, Rcpp::NumericVector r, arma::mat W) {

  double lb[4] = {0, 0, 0, 0};                               /* lower bounds */
  double ub[4] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, HUGE_VAL};   /* upper bounds */
  nlopt_opt opt;
  opt = nlopt_create(NLOPT_LN_NELDERMEAD, 4);                /* algorithm and dimensionality */
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  Rcpp::List lista = Rcpp::List::create(Rcpp::Named("delta") = Delta,
                                        Rcpp::Named("r") = r,
                                        Rcpp::Named("W") = W);
  my_data &temp = (my_data &) lista;
  nlopt_set_min_objective(opt, Q1_w, &temp);

  nlopt_set_xtol_rel(opt, 1e-8);
  double minf;                                                /* the minimum objective value, upon return */
  // nlopt_result result = nlopt_optimize(opt, x.begin(), &minf);
  nlopt_optimize(opt, x.begin(), &minf);
  nlopt_destroy(opt);                                         /* dispose of the nlopt_opt object */

  Rcpp::NumericVector theta(4);
  for (int i = 0; i < 4; ++i) theta[i] = x[i];

  double res = minf;
  // Rcpp::Rcout << "nlopt_result: " << result << "\n";

  return Rcpp::List::create(Rcpp::Named("par") = theta, Rcpp::Named("val") = res);
}


// [[Rcpp::export]]
Rcpp::List GMM_CKLS(Rcpp::NumericVector r, double delta, Rcpp::NumericVector guess, int maxiter, double tol1, double tol2) {

  int n = r.size(), iter = 0;
  int ell = n - 2;
  double res;
  Rcpp::NumericVector w(ell), theta_1(4), theta_1_tmp(4), theta_2(4);
  arma::mat W, Shat;
  Rcpp::List fit2;

  IntegerVector tmp = Rcpp::Range(1, ell);
  w = 1 - as<Rcpp::NumericVector>(tmp) / (ell + 1);
  Rcpp::List fit1 = optQ(guess, delta, r);
  theta_1 = fit1[0];
  res = fit1[1];

  bool seguir = true;
  while (seguir) {
    iter += 1;
    // Rcpp::Rcout << "iter: " << iter << "\n";

    Shat = S(theta_1[0], theta_1[1], theta_1[2], theta_1[3], delta, r, 0, 0); // 4x4
    for (int i = 0; i < ell; ++i) {
      Shat += w[i] * (S(theta_1[0], theta_1[1], theta_1[2], theta_1[3], delta, r, i+1, i+1) +
              trans(S(theta_1[0], theta_1[1], theta_1[2], theta_1[3], delta, r, i+1, i+1)));
    }

    W = Shat.i();

    theta_1_tmp = Rcpp::clone(theta_1);
    fit2        = optQ1(theta_1, delta, r, W);  // Pass by reference, changes the value of theta_1
    theta_2     = fit2[0];
    res         = fit2[1];

    if ((sum(abs(theta_1_tmp - theta_2)) < tol1) | (res < tol2) | (iter > maxiter)) {
      seguir = false;
    }
    theta_1 = theta_2;
    // Rcpp::Rcout << "theta1: " << theta_1 << "\n\n";
  }
  return Rcpp::List::create(Rcpp::Named("par") = theta_1, Rcpp::Named("W") = W);
}


