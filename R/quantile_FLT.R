#' Quantile Functional Linear Test (FLT) for GLS Spectrum Regression
#'
#' Performs a Wald-type integral test to assess whether the estimated coefficient curve
#' \eqn{\beta_1(t)} is significantly different from zero over a specified quantile interval
#' \eqn{[a, b]}. This test is based on the output of \code{\link{wasserstein_spectrum}}.
#'
#' @param spectrum_obj An object returned by \code{\link{wasserstein_spectrum}}, which contains
#' basis matrix, coefficient estimates, robust covariance matrix, and other fitting information.
#' @param a Lower bound of the quantile interval (e.g., 0.25).
#' @param b Upper bound of the quantile interval (e.g., 0.75).
#'
#' @return A list containing:
#' \item{a}{Lower bound of the tested quantile interval.}
#' \item{b}{Upper bound of the tested quantile interval.}
#' \item{delta}{Estimated integral of \(\hat{\beta}_1(t)\) over \([a, b]\).}
#' \item{se}{Estimated standard error of the integral.}
#' \item{z}{Wald test statistic.}
#' \item{p}{Two-sided p-value for the test.}
#'
#' @details
#' The function computes the test statistic:
#' \deqn{
#' Z = \frac{\int_a^b \hat{\beta}_1(t) dt}{\sqrt{\mathrm{Var}\left( \int_a^b \hat{\beta}_1(t) dt \right)}}
#' }
#' using the B-spline basis expansion and robust sandwich variance estimates from the
#' generalized least squares spectrum fit.
#'
#' The integral is approximated using a Riemann sum over the B-spline basis,
#' and the variance of the integral is computed using the delta method.
#'
#' @examples
#' \dontrun{
#' fit <- wasserstein_spectrum(df = example_df,
#'                                 feature_col = "Shannon",
#'                                 outcome_col = "BMI",
#'                                 confounder_cols = c("Age", "Sex"))
#' flt_result <- quantile_FLT(fit, a = 0.25, b = 0.75)
#' print(flt_result$p)
#' }
#'
#' @importFrom stats pnorm
#' @export
quantile_FLT <- function(spectrum_obj, a, b) {
  # Extract necessary info
  Phi <- spectrum_obj$basis          # m x K
  t_grid <- spectrum_obj$quantiles   # length m
  coef <- spectrum_obj$coef          # length p*K
  vcov <- spectrum_obj$vcov          # (p*K) x (p*K)
  idx_x <- spectrum_obj$idx_x        # scalar: column index of outcome var
  p <- spectrum_obj$basis_df         # spline df = K
  X_cols <- spectrum_obj$idx_x
  X0_dim <- length(coef) / p
  
  # Subset quantiles in [a, b]
  idx_ab <- which(t_grid >= a & t_grid <= b)
  Phi_sub <- Phi[idx_ab, , drop = FALSE]  # subgrid
  delta_t <- mean(diff(t_grid))           # assume uniform grid
  
  # Compute c vector: ∫_a^b phi(τ) dτ ≈ sum φ(τ_j) * Δτ
  c_vec <- colSums(Phi_sub) * delta_t     # length K
  
  # Select alpha coefficients: positions of outcome predictor across basis
  sel_idx <- idx_x + (0:(p - 1)) * X0_dim
  alpha_hat <- coef[sel_idx]
  cov_alpha <- vcov[sel_idx, sel_idx]
  
  # Test statistic
  delta_hat <- sum(c_vec * alpha_hat)
  var_delta <- as.numeric(t(c_vec) %*% cov_alpha %*% c_vec)
  z_score <- delta_hat / sqrt(var_delta)
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  return(list(
    a = a,
    b = b,
    delta = delta_hat,
    se = sqrt(var_delta),
    z = z_score,
    p = p_value
  ))
}
