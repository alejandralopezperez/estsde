#include <Rcpp.h>

using namespace Rcpp;

Rcpp::List f_OU(int , Rcpp::NumericVector & , double , Rcpp::NumericVector , int );
void       LDLt(int , Rcpp::NumericMatrix & );
double     loglike(int , Rcpp::NumericVector , double , Rcpp::NumericVector );


//[[Rcpp::export]]
Rcpp::List MCMC_ou(const Rcpp::NumericVector& obs, int n_obs, int m_data, int niter,
                   int total_iter, double alpha, double kappa, double sigma) {

  int                  M = 3, j, k, bbl, nst, iter;

  double               DT = 1.0, dt, tmp, tmp1, p1, p2,
                       p3, sum, logprob, xold, xnew,
                       dx, L00, L10, L11, es1, es2, es3;

  Rcpp::NumericVector  phat(M),
                       link(m_data+1),
                       randn(m_data+1),
                       mean(m_data+1),
                       y(n_obs*m_data+1);

  Rcpp::NumericMatrix  param(total_iter, M),
                       chain(n_obs, m_data+1),
                       L(m_data-1, m_data-1),
                       X(n_obs*m_data+1, 2),
                       XXTinv(2,2);

  param(0,0) = alpha;
  param(0,1) = kappa;
  param(0,2) = sigma;

  dt = DT/m_data;
  xnew = obs(0);
  for (k=0; k<n_obs; k++ ) {
    chain(k,0) = xnew;
    xold = xnew;
    xnew = obs(k+1);
    dx = (xnew-xold)/m_data;
    for ( j=1 ; j<m_data ; j++ ) {
      chain(k,j) = xold+dx*j;
    }
    chain(k,m_data) = xnew;
  }

  iter = 1;
  while ( iter < total_iter ) {

    for (k=0; k<M; k++) phat[k] = param(iter-1,k);

    for (k=0; k<n_obs; k++) {

      bbl=0;

      while (bbl < m_data-1) {

        do {
          nst = R::rpois(1.0)+1;
        } while (nst <= 1);

        if (bbl+nst > m_data) nst = m_data-bbl;
        link(0) = mean(0) = chain(k,bbl);
        link(nst) = mean(nst) = chain(k,bbl+nst);
        dx = (mean(nst)-mean(0))/((double) nst);
        for (j=1; j<nst; j++) {
          mean(j) = link(0)+dx*((double) j);
        }

        Rcpp::List rdo = f_OU(nst-1, mean, dt, phat, m_data);
        Rcpp::NumericMatrix V = rdo(0);
        Rcpp::NumericVector u = rdo(1);

        LDLt(nst-1, V);
        for (j=0; j<nst-2; j++) {
          tmp = sqrt( std::abs(V(j,j)) );
          L(j,j) = tmp;
          L(j+1,j) = V(j+1,j)*tmp;
        }
        L(nst-2,nst-2) = sqrt( std::abs(V(nst-2,nst-2)) );

        for (j=1; j<nst; j++) {
          randn(j) = R::rnorm(0.0,1.0);
        }

        link(nst-1) = randn(nst-1) / L(nst-2,nst-2);
        for (j=nst-2; j>0; j--) {
          link(j) = (randn(j)-L(j,j-1)*link(j+1)) / L(j-1,j-1);
        }

        for (j=1; j<nst; j++) {
          link(j) += mean(j);
        }

        // new
        logprob = -loglike(nst, link, dt, phat);

        // old
        u(0) = link(0);
        u(nst) = link(nst);
        for (j=1; j<nst; j++) {
          u(j) = chain(k,bbl+j);
        }
        logprob += loglike(nst, u, dt, phat);

        sum = 0.0;
        for (j=1; j<nst; j++) {
          sum += randn(j)*randn(j);
        }
        logprob += 0.5*sum;

        for (j=1; j<nst; j++) {
          u(j-1) = chain(k,bbl+j)-mean(j);
        }

        for (j=0; j<nst-2; j++) {
          randn(j+1) = L(j,j)*u(j)+L(j+1,j)*u(j+1);
        }
        randn(nst-1) = L(nst-2,nst-2)*u(nst-2);

        sum = 0.0;
        for (j=1; j<nst; j++) {
          sum += randn(j)*randn(j);
        }
        logprob -= 0.5*sum;

        // change
        if ( logprob >= 0.0 ) {
          for (j=1; j<nst; j++) {
            chain(k,bbl+j) = link(j);
          }
        } else {
          tmp = std::exp(logprob);
          if ( R::runif(0,1) < tmp ) {
            for (j=1; j<nst; j++) {
              chain(k,bbl+j) = link(j);
            }
          }
        }
        bbl += nst-1;
      }
    }


    for (k=0; k<n_obs; k++) {
      for (j=0; j<m_data; j++) {
        X(k*m_data+j,0) =  pow(dt,0.5);
        X(k*m_data+j,1) = -pow(dt,0.5) * pow(chain(k,j), 1.0);
        y(k*m_data+j) = (chain(k,j+1)-chain(k,j)) / pow(dt,0.5);
      }
    }

    XXTinv(0,0) = XXTinv(0,1) = XXTinv(1,1) = 0.0;
    for (j=0; j<n_obs*m_data; j++) {
      XXTinv(1,1) += X(j,0)*X(j,0);
      XXTinv(0,1) -= X(j,0)*X(j,1);
      XXTinv(0,0) += X(j,1)*X(j,1);
    }

    tmp = 1.0/(XXTinv(0,0)*XXTinv(1,1)-XXTinv(0,1)*XXTinv(0,1));
    XXTinv(0,0) *= tmp;
    XXTinv(1,1) *= tmp;
    XXTinv(0,1) *= tmp;
    XXTinv(1,0) = XXTinv(0,1);

    for (j=0; j<2; j++) {
      randn(j) = R::rnorm(0.0,1.0);
    }

    tmp = 0.0;
    tmp1 = 0.0;
    for (j=0 ; j<n_obs*m_data ; j++ ) {
      tmp  += X(j,0)*y(j);
      tmp1 += X(j,1)*y(j);
    }
    mean[0] = XXTinv(0,0)*tmp + XXTinv(0,1)*tmp1;
    mean[1] = XXTinv(1,0)*tmp + XXTinv(1,1)*tmp1;

    tmp = 0.0;
    for (j=0; j<n_obs*m_data; j++) {
      tmp1 = y(j) - X(j,0)*mean(0) - X(j,1)*mean(1);
      tmp += tmp1*tmp1;
    }

    param(iter,2) = sqrt( 1.0/R::rgamma(0.5*((double) n_obs*m_data-2), 2.0/tmp) );

    tmp = pow(XXTinv(0,0),0.5);
    L00 = param(iter,2) * tmp;
    L10 = param(iter,2) * XXTinv(1,0)/tmp;
    L11 = param(iter,2) * pow(XXTinv(0,0)*XXTinv(1,1)-XXTinv(0,1)*XXTinv(0,1), 0.5) / tmp;

    param(iter,0) = mean(0) + L00*randn(0);
    param(iter,1) = mean(1) + L10*randn(0) + L11*randn(1);

    iter++;
  }

  //-------------------------------------------------
  p1 = p2 = p3 = 0.0;
  for ( k=(total_iter-niter) ; k<total_iter ; k++ ) {
    p1 += param(k,0);
    p2 += param(k,1);
    p3 += param(k,2);
  }
  // Media
  p1 /= ((double) niter);
  p2 /= ((double) niter);
  p3 /= ((double) niter);

  es1 = es2 = es3 = 0.0;
  for ( k=(total_iter-niter) ; k<total_iter ; k++ ) {
    es1 += (param(k,0)-p1) * (param(k,0)-p1);
    es2 += (param(k,1)-p2) * (param(k,1)-p2);
    es3 += (param(k,2)-p3) * (param(k,2)-p3);
  }
  es1 /= ((double) niter-1);
  es2 /= ((double) niter-1);
  es3 /= ((double) niter-1);

  // Error estandar
  es1 = sqrt(es1);
  es2 = sqrt(es2);
  es3 = sqrt(es3);

  Rcpp::NumericVector theta = Rcpp::NumericVector::create(p1, p2, p3);
  Rcpp::NumericVector error = Rcpp::NumericVector::create(es1, es2, es3);

  return Rcpp::List::create(Rcpp::Named("param") = theta,
                                      _("error") = error);
}


Rcpp::List f_OU(int len, Rcpp::NumericVector &mean, double dt, Rcpp::NumericVector phat, int m_data) {
  int i, j;
  double c, g, gold, d, dold;
  Rcpp::NumericVector u(m_data+1), tmpmean(len), lsub(len-1), umain(len), usup(len-1);
  Rcpp::NumericMatrix V(m_data-1, m_data-1);

  c = 1.0-phat(1)*dt;

  for (i=0; i<5; i++) {

    std::fill(V.begin(), V.end(), 0.0);

    gold = 1.0/(phat(2)*phat(2));
    dold = mean(1)-mean(0)-(phat(0)-phat(1)*mean(0));
    for (j=1; j<=len; j++) {

      g = 1.0/(phat(2)*phat(2));
      d = mean(j+1)-mean(j)-(phat(0)-phat(1)*mean(j));

      u(j-1) = -1.0/dt * (gold*dold - g*d*c);
      V(j-1,j-1) = 1.0/dt * (gold + g*c*c);
      if (j == len) {
      } else {
        V(j-1,j) = V(j,j-1) = 1.0/dt * (-g*c);  // Fuera de la diagonal
      }

      gold = g;
      dold = d;
    }


    for (j=0; j<len-1; j++) {
      usup(j) = V(j,j+1);
    }
    umain(0) = V(0,0);
    for (j=1; j<len; j++) {
      lsub(j-1) = V(j,j-1)/umain(j-1);
      umain(j)  = V(j,j)-lsub(j-1)*usup(j-1);
    }
    tmpmean(0) = u(0);
    for (j=1; j<len; j++) {
      tmpmean(j) = u(j)-lsub(j-1)*tmpmean(j-1);
    }
    u(len-1) = tmpmean(len-1)/umain(len-1);
    for (j=len-2; j>=0; j--) {
      u(j) = (tmpmean(j)-usup(j)*u(j+1))/umain(j);
    }

    for (j=0; j<len; j++) {
      mean(j+1) += u(j);
    }
  }

  gold = 1.0/(phat[2]*phat[2]);

  for (j=1; j<=len; j++) {

    g = 1.0/(phat[2]*phat[2]);

    V(j-1,j-1) = 1.0/dt * (gold + g*c*c);
    if (j == len) {
    } else {
      V(j-1,j) = V(j,j-1) = 1.0/dt * (-g*c);
    }
  }

  return Rcpp::List::create(Rcpp::Named("V") = V,
                            _("u") = u);
}


void LDLt(int n, Rcpp::NumericMatrix &a) {
  double  sum, pivot=0.0;
  int     i, j, k;

  if ( a(0,0) != 0.0 ) {
    pivot = 1.0/a(0,0);
  } else {
    Rcpp::Rcout << "\nSingular matrix"<< "\n";
  }
  for (k=1 ; k<n ; k++) a(k,0) *= pivot;

  for (i=1 ; i<n ; i++) {
    for (sum=0.0, k=0; k<i; k++) sum += a(k,k)*a(i,k)*a(i,k);
    a(i,i) -= sum;
    if ( a(i,i) != 0.0 ) {
      pivot = 1.0/a(i,i);
    } else {
      Rcpp::Rcout << "\nSingular matrix"<< "\n";
    }
    for ( j=i+1 ; j<n ; j++ ) {
      for ( sum=0.0, k=0 ; k<i ; k++ ) sum += a(k,k)*a(i,k)*a(j,k);
      a(j,i) = (a(j,i)-sum)*pivot;
    }
  }
  return;
}


double loglike(int nst, Rcpp::NumericVector link, double dt, Rcpp::NumericVector phat) {
  int    i;
  double kappa, alpha, sigma, tmp, negll, pi = M_PI;

  alpha = phat(0);
  kappa = phat(1);
  sigma = phat(2);

  negll = 0.5*((double) nst)*log(2.0*sigma*sigma*pi*dt);
  for (i=1; i <= nst; i++) {
    tmp = link[i]-link[i-1]-(alpha-kappa*link[i-1])*dt;
    negll += 0.5*(tmp*tmp/(sigma*sigma*dt));
  }
  return negll;
}

