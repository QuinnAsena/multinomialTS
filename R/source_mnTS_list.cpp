#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List mnTS_ml_cpp_list(
    arma::vec par,
    arma::vec par_fixed,
    arma::mat Y,
    arma::mat X,
    arma::vec Tsample) {

    Tsample -= 1;

    arma::vec parameters = par_fixed;

    arma::uvec par_nas = arma::find_nonfinite(par_fixed);
    parameters.elem(par_nas) = par;

    int n = Y.n_rows;

    size_t Tmax = static_cast<size_t>(Tsample.back());

    arma::mat B0(n, 1);
    B0.col(0) = parameters.head(n);

    arma::vec subvector_C = parameters.subvec(n, n + n * n - 1);
    arma::mat C = arma::reshape(subvector_C, n, n);

    double sigma = parameters(n + n * n);

    arma::vec subvector_L = parameters.subvec(n + n * n + 1, n + n * n + n * n);
    arma::mat L = arma::reshape(subvector_L, n, n);
    arma::mat S = arma::trans(L) * L;

    double Ssd = S(0, 0);

    double sigmaSquared = sigma * sigma;
    arma::mat Se = sigmaSquared * S / Ssd;

    double dispersion = std::exp(parameters(n + n * n + 1 + n * n));

    arma::mat B;
    if (X.n_rows > 0) {
      arma::uword p = X.n_rows;
      arma::vec subvector_B = parameters.subvec(n + n * n + 1 + n * n + 1, n + n * n + 1 + n * n + p * n);
      B = arma::trans(arma::reshape(subvector_B, p, n));
    }

    double logFt = 0;
    double vFv = 0;

    arma::mat y;
    arma::mat PP;

    arma::mat y_fitted, mu_fitted, se_y_fitted, se_mu_upper_fitted, se_mu_lower_fitted;

    y_fitted = arma::mat(n, Tmax+1, arma::fill::none);  // Initialize with NaN values
    mu_fitted = arma::mat(n, Tmax+1, arma::fill::none);
    // std::cout << "mu_fitted =\n " << mu_fitted << std::endl;

    se_y_fitted = arma::mat(n, Tmax+1, arma::fill::none);
    se_mu_upper_fitted = arma::mat(n, Tmax+1, arma::fill::none);
    se_mu_lower_fitted = arma::mat(n, Tmax+1, arma::fill::none);


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

      arma::vec mu;
      if (!Y.cols(arma::find(t == Tsample)).has_nan() & arma::any(t == Tsample)) {
        double size = arma::accu(Y.cols(arma::find(t == Tsample)));
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
          double LL = 1e+10;
          return Rcpp::List::create(
            Rcpp::Named("LL") = LL
          );
          }

        arma::mat invF = arma::inv(FF);
        arma::mat v = Y(arma::regspace<arma::uvec>(1, Y.n_rows - 1), arma::find(t == Tsample));
        v /= size;
        v -= mu.subvec(1, mu.n_elem -1);

        y += PP * ZDM.t() * invF * v;

        PP -= PP * ZDM.t() * invF * ZDM * PP;

        logFt += std::log(std::abs(arma::det(FF)));
        vFv += arma::as_scalar(v.t() * invF * v);
      }

      // Update 'y.fitted' column 't' with 'y'
      y_fitted.col(t) = y;

      // Update 'mu.fitted' column 't' with 'mu'
      mu = arma::exp(y);
      mu /= arma::sum(mu);

      mu_fitted.col(t) = mu;

      // Calculate 'se.y.fitted' as the square root of the diagonal of 'PP'
      se_y_fitted.col(t) = arma::sqrt(arma::diagvec(PP));

      se_mu_upper_fitted.col(t) = arma::exp(y_fitted.col(t) + se_y_fitted.col(t))/arma::as_scalar(arma::sum(arma::exp(y)));
      se_mu_lower_fitted.col(t) = arma::exp(y_fitted.col(t) - se_y_fitted.col(t))/arma::as_scalar(arma::sum(arma::exp(y)));

  }
    double LL = logFt + vFv;

    if (!arma::is_finite(LL)) {
      LL = 1e+10;

      return Rcpp::List::create(
        Rcpp::Named("LL") = LL
      );

    } else {
      arma::mat V = arma::trans(L) * L;
      arma::mat par = parameters.elem(arma::find_nan(par_fixed));

      int Tsample_len = Tsample.n_elem;

      double logLik = -((n - 1) * Tsample_len / 2) * std::log(2 * arma::datum::pi) - LL / 2;

      return Rcpp::List::create(
        Rcpp::Named("LL") = LL /= 2,
        Rcpp::Named("B0") = B0,
        Rcpp::Named("C") = C,
        Rcpp::Named("sigma") = std::abs(sigma),
        Rcpp::Named("V") = V,
        Rcpp::Named("B") = B,
        Rcpp::Named("dispersion") = dispersion,
        Rcpp::Named("logLik") = logLik,
        Rcpp::Named("par") = par,
        Rcpp::Named("mu") = arma::trans(mu_fitted),
        Rcpp::Named("y") = arma::trans(y_fitted),
        Rcpp::Named("se.y.fitted") = arma::trans(se_y_fitted),
        Rcpp::Named("se.mu.upper.fitted") = arma::trans(se_mu_upper_fitted),
        Rcpp::Named("se.mu.lower.fitted") = arma::trans(se_mu_lower_fitted)
      );
    }
  }
