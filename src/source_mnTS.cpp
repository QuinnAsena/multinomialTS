#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double mnTS_ml_cpp(
    arma::vec par,
    arma::vec par_fixed,
    arma::mat Y,
    arma::mat X,
    arma::vec Tsample) {

    Tsample -= 1;
  // std::cout << "Tsample =\n " << Tsample << std::endl;
    arma::vec parameters = par_fixed;

    arma::uvec par_nas = arma::find_nonfinite(par_fixed);
    parameters.elem(par_nas) = par;

    arma::uword n = Y.n_rows;

    size_t Tmax = static_cast<size_t>(Tsample.back());
    // std::cout << "Tmax =\n " << Tmax << std::endl;
    arma::mat B0(n, 1);
    B0.col(0) = parameters.head(n);

    arma::vec subvector_C = parameters.subvec(n, n + n * n - 1);
    arma::mat C = arma::reshape(subvector_C, n, n);

    double sigma = parameters(n + n * n);
    // std::cout << "sigma =\n " << sigma << std::endl;

    arma::vec subvector_L = parameters.subvec(n + n * n + 1, n + n * n + n * n);
    arma::mat L = arma::reshape(subvector_L, n, n);
    arma::mat S = arma::trans(L) * L;

    double Ssd = S(0, 0);

    double sigmaSquared = sigma * sigma;
    arma::mat Se = sigmaSquared * S / Ssd;
    // std::cout << "Se =\n " << Se << std::endl;
    // std::cout << "L =\n " << L << std::endl;


    double dispersion = std::exp(parameters(n + n * n + 1 + n * n));

    arma::mat B;
    if (X.n_rows > 0) {
      arma::uword p = X.n_rows;
      arma::vec subvector_B = parameters.subvec(n + n * n + 1 + n * n + 1, n + n * n + 1 + n * n + p * n);
      // std::cout << "subvector_B =\n " << subvector_B << std::endl;
      // std::cout << "subvector_B =\n " << subvector_B << std::endl;
      B = arma::trans(arma::reshape(subvector_B, p, n));
    }
    // std::cout << "B0 =\n " << B0 << std::endl;
    // std::cout << "B =\n " << B << std::endl;
    // std::cout << "C =\n " << C << std::endl;
    // std::cout << "sigma =\n " << sigma << std::endl;
    // std::cout << "L =\n " << L << std::endl;
    // std::cout << "S =\n " << S << std::endl;

    double logFt = 0;
    double vFv = 0;

    arma::mat y;
    arma::mat PP;
    // CHECK INDEXING edited to Tmax + 1
    // CHECK INDEXING Check against Tsample
    // CHECK INDEXING Seems to be working but:
    // CHECK INDEXING!!!!

    // if (arma::max(arma::abs(arma::eig_gen(C))) > 1) {
    //   double LL = 1e10;
    //   return LL;
    // }

    for (arma::uword t = 0; t < Tmax + 1; ++t) {
      if (t == 0) {
        // Initial values at the stationary distribution
        if (!X.is_empty()) {
          y = B0 + B * X.col(t);
        } else {
          y = B0;
        }

        if (arma::max(arma::abs(arma::eig_gen(C))) < 1) {
          arma::mat Kronecker = arma::kron(C, C);
          arma::mat I = arma::eye(n * n, n * n);
          PP = arma::solve(I - Kronecker, arma::vectorise(Se));
          PP = arma::reshape(PP, n, n);
        } else {
          PP = Se;
        }
      } else {
        // PREDICTION EQUATIONS
        if (!X.is_empty()) {
          y = B0 + C * (y - B0 - B * X.col(t - 1)) + B * X.col(t);
        } else {
          y = B0 + C * (y - B0);
        }
        PP = C * PP * C.t() + Se;
      }

      // std::cout << "Tsample =\n " << Tsample << std::endl;
      arma::vec mu;
      if (!Y.cols(arma::find(t == Tsample)).has_nan() & arma::any(t == Tsample)) {
        // std::cout << "t =\n " << t << std::endl;
        double size = arma::accu(Y.cols(arma::find(t == Tsample)));
        // std::cout << "size =\n " << size << std::endl;
        mu = arma::exp(y);
        mu /= arma::sum(mu);

        arma::mat W = arma::diagmat(arma::pow(mu % (1 - mu), -1)) * (arma::diagmat(mu) - mu * mu.t()) * arma::diagmat(arma::pow(mu % (1 - mu), -1));
        arma::mat M = arma::diagmat(1 + mu) - (arma::ones<arma::mat>(n, 1) * mu.t());
        arma::mat D = arma::diagmat(mu % (1 - mu));
        arma::mat ZDM = (D * M);
        ZDM = ZDM.rows(1, ZDM.n_rows - 1);

        arma::mat intermediate = D * (W / size) * D;
        arma::mat submatrix = intermediate.submat(1, 1, intermediate.n_rows - 1, intermediate.n_cols - 1);
        arma::mat FF = ZDM * PP * ZDM.t() + dispersion * submatrix;

        if (FF.has_nan() || arma::rcond(FF) < 1e-12) {
          return 1e+10;
          }

        arma::mat invF = arma::inv(FF);
        arma::mat v = Y(arma::regspace<arma::uvec>(1, Y.n_rows - 1), arma::find(t == Tsample));
        v /= size;
        v -= mu.subvec(1, mu.n_elem -1);

        y += PP * ZDM.t() * invF * v;
        PP -= PP * ZDM.t() * invF * ZDM * PP;

        logFt += std::log(std::abs(arma::det(FF)));
        vFv += arma::as_scalar(v.t() * invF * v);
        // std::cout << "logFt =\n " << logFt << std::endl;

        // arma::vec vector = arma::linspace<arma::vec>(129, 138, 10);
        // if (arma::any(t == vector)) {
          // std::cout << "t =\n " << t << std::endl;
        //   std::cout << "Tsample =\n " << arma::find(t == Tsample) << std::endl;
        //   std::cout << "logFt =\n " << logFt << std::endl;
        //   std::cout << "vFv =\n " << vFv << std::endl;
        // }

      }
  }
    double LL = logFt + vFv;

    if (!arma::is_finite(LL)) {
      LL = 1e+10;

      return LL;

    } else {
      // std::cout << "LL =\n " << LL/2 << std::endl;

      return LL /= 2;
    }
  }
