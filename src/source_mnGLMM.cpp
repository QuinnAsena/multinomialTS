#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List mnGLMM_ml_cpp_list(
    arma::vec par,
    arma::vec par_fixed,
    arma::mat Y,
    arma::mat X) {

    // Free parameters are the NaN entries of par_fixed. This must match R's
    // is.na(par.fixed), which is NaN only: find_nonfinite would also pick up
    // +/-Inf entries and silently shift the parameter mapping.
    arma::vec parameters = par_fixed;
    arma::uvec par_nas = arma::find_nan(par_fixed);
    parameters.elem(par_nas) = par;

    arma::uword n    = Y.n_cols;
    arma::uword p    = X.n_cols;
    arma::uword Tmax = Y.n_rows;

    // --- Parameter extraction ---
    arma::mat B = arma::reshape(parameters.head(n * p), p, n);

    arma::mat L   = arma::reshape(parameters.subvec(n * p, n * p + n * n - 1), n, n);
    arma::mat V   = arma::trans(L) * L;
    double    Ssd = V(0, 0);
    double    sigma = parameters(n * p + n * n);
    arma::mat P   = (sigma * sigma) * V / Ssd;

    double dispersion = std::exp(parameters(n * p + n * n + 1));

    arma::mat y_mat    = X * B;
    arma::mat exp_y    = arma::exp(y_mat);
    arma::vec row_sums = arma::sum(exp_y, 1);
    arma::mat mu_mat   = exp_y.each_col() / row_sums;

    arma::vec size_vec = arma::sum(Y, 1);

    double SS = 0.0;

    for (arma::uword i = 0; i < Tmax; ++i) {

        arma::vec mu_i  = mu_mat.row(i).t();
        double   size_i = size_vec(i);

        arma::vec w_diag = arma::pow(mu_i % (1.0 - mu_i), -1.0);
        arma::mat W = arma::diagmat(w_diag)
                    * (arma::diagmat(mu_i) - mu_i * mu_i.t())
                    * arma::diagmat(w_diag);

        arma::mat M     = arma::diagmat(1.0 + mu_i)
                        - arma::ones<arma::mat>(n, 1) * mu_i.t();
        arma::mat S_cov = dispersion * (W / size_i) + M * P * M.t();

        arma::mat D = arma::diagmat(mu_i % (1.0 - mu_i));
        S_cov = D * S_cov * D;

        arma::vec v = arma::vectorise(Y.submat(i, 1, i, n - 1)) / size_i
                    - mu_i.subvec(1, n - 1);

        arma::mat S_sub = S_cov.submat(1, 1, n - 1, n - 1);

        if (S_sub.has_nan() || arma::rcond(S_sub) < 1e-12)
            return Rcpp::List::create(Rcpp::Named("LL") = 1e+10);

        arma::mat invS = arma::inv(S_sub);

        double log_det_val, sign;
        arma::log_det(log_det_val, sign, S_sub);

        SS += arma::as_scalar(v.t() * invS * v) + log_det_val;
    }

    if (!std::isfinite(SS))
        return Rcpp::List::create(Rcpp::Named("LL") = 1e+10);

    arma::vec par_out = parameters.elem(arma::find_nan(par_fixed));
    double logLik = -((double)(n - 1) * Tmax / 2.0) * std::log(2.0 * arma::datum::pi) - SS / 2.0;

    return Rcpp::List::create(
        Rcpp::Named("LL")         = SS / 2.0,
        Rcpp::Named("B")          = B,
        Rcpp::Named("sigma")      = std::abs(sigma),
        Rcpp::Named("V")          = V,
        Rcpp::Named("dispersion") = dispersion,
        Rcpp::Named("logLik")     = logLik,
        Rcpp::Named("par")        = par_out,
        Rcpp::Named("mu")         = mu_mat
    );
}
