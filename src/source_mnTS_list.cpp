#include <RcppArmadillo.h>
#include <Rcpp.h>

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

    // Convert Tsample to 0-indexed to match C++ loop
    Tsample -= 1;

    // Free parameters are the NaN entries of par_fixed. This must match R's
    // is.na(par.fixed), which is NaN only: find_nonfinite would also pick up
    // +/-Inf entries and silently shift the parameter mapping.
    arma::vec parameters = par_fixed;
    arma::uvec par_nas = arma::find_nan(par_fixed);
    parameters.elem(par_nas) = par;

    int n = Y.n_rows;
    // max(), not back(): an unsorted Tsample would otherwise truncate the run.
    size_t Tmax = static_cast<size_t>(Tsample.max());

    // --- Parameter extraction ---
    arma::mat B0(n, 1);
    B0.col(0) = parameters.head(n);

    arma::mat C = arma::reshape(parameters.subvec(n, n + n * n - 1), n, n);

    double sigma = parameters(n + n * n);

    arma::mat L = arma::reshape(parameters.subvec(n + n * n + 1, n + n * n + n * n), n, n);
    arma::mat S = arma::trans(L) * L;
    double Ssd = S(0, 0);
    arma::mat Se = (sigma * sigma) * S / Ssd;

    double dispersion = std::exp(parameters(n + n * n + 1 + n * n));

    // Use the same emptiness test as the prediction equations below, so there
    // is one definition of "no covariates" throughout.
    arma::mat B;
    if (!X.is_empty()) {
        arma::uword p = X.n_rows;
        arma::vec subvector_B = parameters.subvec(n + n * n + 1 + n * n + 1,
                                                   n + n * n + 1 + n * n + p * n);
        B = arma::trans(arma::reshape(subvector_B, p, n));
    }

    double logFt = 0.0;
    double vFv   = 0.0;

    arma::mat y;
    arma::mat PP;

    // Initialise fitted-value matrices with NaN to match R's matrix(NA, ...)
    arma::mat y_fitted         (n, Tmax + 1, arma::fill::value(arma::datum::nan));
    arma::mat mu_fitted        (n, Tmax + 1, arma::fill::value(arma::datum::nan));
    arma::mat se_y_fitted      (n, Tmax + 1, arma::fill::value(arma::datum::nan));
    arma::mat se_mu_upper_fitted(n, Tmax + 1, arma::fill::value(arma::datum::nan));
    arma::mat se_mu_lower_fitted(n, Tmax + 1, arma::fill::value(arma::datum::nan));

    for (arma::uword t = 0; t < Tmax + 1; ++t) {

        if (t == 0) {
            // Initial state at the stationary distribution
            y = X.is_empty() ? B0 : B0 + B * X.col(t);

            if (arma::max(arma::abs(arma::eig_gen(C))) < 1) {
                arma::mat Kronecker = arma::kron(C, C);
                arma::mat I = arma::eye(n * n, n * n);
                PP = arma::reshape(arma::solve(I - Kronecker, arma::vectorise(Se)), n, n);
            } else {
                PP = Se;
            }

        } else {
            // Prediction equations
            if (X.is_empty()) {
                y = B0 + C * (y - B0);
            } else {
                y = B0 + C * (y - B0 - B * X.col(t - 1)) + B * X.col(t);
            }
            PP = C * PP * C.t() + Se;
        }

        // Updating equations — only when t is an observation time with no missing data
        arma::uvec t_idx = arma::find(Tsample == static_cast<double>(t));
        if (t_idx.n_elem > 0 && !Y.cols(t_idx).has_nan()) {

            double size = arma::accu(Y.cols(t_idx));

            arma::vec mu = arma::exp(y);
            mu /= arma::sum(mu);

            arma::mat W   = arma::diagmat(arma::pow(mu % (1 - mu), -1))
                          * (arma::diagmat(mu) - mu * mu.t())
                          * arma::diagmat(arma::pow(mu % (1 - mu), -1));
            arma::mat M   = arma::diagmat(1 + mu) - (arma::ones<arma::mat>(n, 1) * mu.t());
            arma::mat D   = arma::diagmat(mu % (1 - mu));
            arma::mat DM  = D * M;
            arma::mat ZDM = DM.rows(1, n - 1);

            arma::mat intermediate = D * (W / size) * D;
            arma::mat FF = ZDM * PP * ZDM.t()
                         + dispersion * intermediate.submat(1, 1, n - 1, n - 1);

            if (FF.has_nan() || arma::rcond(FF) < 1e-12) {
                return Rcpp::List::create(Rcpp::Named("LL") = 1e+10);
            }

            arma::mat invF = arma::inv(FF);
            arma::mat v = Y(arma::regspace<arma::uvec>(1, n - 1), t_idx);
            v /= size;
            v -= mu.subvec(1, n - 1);

            y  += PP * ZDM.t() * invF * v;
            PP -= PP * ZDM.t() * invF * ZDM * PP;

            // Numerically stable log-determinant
            double log_det_val, sign;
            arma::log_det(log_det_val, sign, FF);
            logFt += log_det_val;

            vFv += arma::as_scalar(v.t() * invF * v);
        }

        // Store fitted values for this time step
        y_fitted.col(t) = y;

        arma::vec mu_t = arma::exp(y);
        mu_t /= arma::sum(mu_t);
        mu_fitted.col(t) = mu_t;

        arma::vec se_y = arma::sqrt(arma::diagvec(PP));
        se_y_fitted.col(t) = se_y;

        double sum_exp_y = arma::as_scalar(arma::sum(arma::exp(y)));
        se_mu_upper_fitted.col(t) = arma::exp(y + se_y) / sum_exp_y;
        se_mu_lower_fitted.col(t) = arma::exp(y - se_y) / sum_exp_y;
    }

    double LL = logFt + vFv;

    if (!std::isfinite(LL)) {
        return Rcpp::List::create(Rcpp::Named("LL") = 1e+10);
    }

    arma::mat V   = arma::trans(L) * L;
    arma::vec par_out = parameters.elem(arma::find_nan(par_fixed));
    int Tsample_len = Tsample.n_elem;
    double logLik = -((n - 1) * Tsample_len / 2.0) * std::log(2 * arma::datum::pi) - LL / 2.0;

    return Rcpp::List::create(
        Rcpp::Named("LL")                  = LL / 2.0,
        Rcpp::Named("B0")                  = B0,
        Rcpp::Named("C")                   = C,
        Rcpp::Named("sigma")               = std::abs(sigma),
        Rcpp::Named("V")                   = V,
        Rcpp::Named("B")                   = B,
        Rcpp::Named("dispersion")          = dispersion,
        Rcpp::Named("logLik")              = logLik,
        Rcpp::Named("par")                 = par_out,
        Rcpp::Named("mu")                  = arma::trans(mu_fitted),
        Rcpp::Named("y")                   = arma::trans(y_fitted),
        Rcpp::Named("se.y.fitted")         = arma::trans(se_y_fitted),
        Rcpp::Named("se.mu.upper.fitted")  = arma::trans(se_mu_upper_fitted),
        Rcpp::Named("se.mu.lower.fitted")  = arma::trans(se_mu_lower_fitted)
    );
}
