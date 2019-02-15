#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute COM-Poisson normalizing constante in log scale
//' @description Computes the normalizing constant in the log scale
//'   usign the LogSumExp trick to avoid numerical issues (See details).
//' @details \code{logspace_add(a, b) = log(exp(a) + exp(b)) = b + log(1 +
//'   exp(a - b))}, where \code{b > a}.
//' @param loglambda A vector of logarithm of the \eqn{\lambda}
//'   parameter.
//' @param nu A vector of dispersion parameters \eqn{\nu}.
//'
//' @return A vector with the logarithms of the resulting
//'   constants.
//' @references Wikipedia, LogSumExp. \url{http://rstudio.com}.
//' @author Eduardo Jr <edujrrib@gmail.com>
//' @export
// [[Rcpp::export]]
NumericVector compute_logz(NumericVector loglambda,
                           NumericVector nu) {
  // Control loop
  int maxiter = 1e4;
  double logepsilon = log(1e-4);
  // Output vector
  int n = nu.size();
  NumericVector out(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logz = 0;
    for (int j = 1; j < maxiter; ++j) {
      double logz_ = j * loglambda[i] - nu[i] * R::lgammafn(j + 1);
      logz = R::logspace_add(logz, logz_);
      if (logz_ < logepsilon) break;
    }
    out[i] = logz;
  }
  return out;
}
