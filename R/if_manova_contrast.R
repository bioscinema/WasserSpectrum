#' Interval-Wise MANOVA Contrast Test for Spectrum Curves
#'
#' Performs a Wald-type contrast test over a specified quantile interval \eqn{[a, b]} to assess whether
#' the integrated group-level spectrum curves differ significantly. This function is based on the output
#' of \code{\link{wasserstein_spectrum_multiclass}}.
#'
#' @param spectrum_obj An object returned by \code{\link{wasserstein_spectrum_multiclass}}, containing
#' estimated group-level coefficients, basis matrix, and robust variance-covariance matrix.
#' @param a Lower bound of the quantile interval \eqn{a \in (0, 1)}.
#' @param b Upper bound of the quantile interval \eqn{b \in (0, 1)}, where \code{b > a}.
#' @param contrast Optional contrast vector or matrix (rows as group comparisons).
#' If \code{NULL} (default), performs a global test across all non-reference groups.
#'
#' @return A list with the following components:
#' \item{statistic}{The Wald statistic (Hotelling-type \eqn{T^2}) for the contrast.}
#' \item{df}{Degrees of freedom for the test (number of contrast rows).}
#' \item{p_value}{P-value from the chi-squared distribution.}
#' \item{contrast}{The contrast matrix used for the test.}
#' \item{estimate}{The estimated contrast effect (linear combination of integrated group curves).}
#' \item{cov_estimate}{Covariance matrix of the estimated contrast effect.}
#'
#' @details
#' This function tests whether the **integrated spectrum curves** differ between groups
#' over the interval \eqn{[a, b]}. The contrast matrix defines the comparison structure
#' (e.g., pairwise or global). If omitted, the function returns a global test comparing
#' all groups versus the reference.
#'
#' The Wald test statistic is:
#' \deqn{
#' T^2 = \theta^\top \mathrm{Var}^{-1}(\theta) \theta,
#' }
#' where \(\theta = C \int_a^b \beta(t)\,dt\).
#'
#' @examples
#' \dontrun{
#' # Global test over middle quantiles
#' fit <- wasserstein_spectrum_multiclass(...)
#' test <- if_manova_contrast(fit, a = 0.2, b = 0.8)
#' 
#' # Pairwise contrast
#' pairwise <- c(1, -1, 0)  # group 1 vs group 2
#' test2 <- if_manova_contrast(fit, a = 0.2, b = 0.8, contrast = pairwise)
#' }
#'
#' @importFrom stats pchisq
#' @export
if_manova_contrast <- function(spectrum_obj, a, b, contrast = NULL) {
  # Extract info
  Phi <- spectrum_obj$basis
  t_grid <- spectrum_obj$quantiles
  coef <- spectrum_obj$coef
  vcov <- spectrum_obj$vcov
  comparison_levels <- spectrum_obj$comparison_levels
  K_group <- length(comparison_levels) + 1
  K_basis <- ncol(Phi)
  p_total <- K_group - 1
  X0_dim <- length(coef) / K_basis
  
  # Compute integration vector c
  idx_ab <- which(t_grid >= a & t_grid <= b)
  Phi_sub <- Phi[idx_ab, , drop = FALSE]
  delta_t <- mean(diff(t_grid))
  c_vec <- colSums(Phi_sub) * delta_t
  
  # Matrix to obtain Delta from alpha
  C_integral <- matrix(0, nrow = p_total, ncol = p_total * K_basis)
  for (k in 1:p_total) {
    start_col <- (k - 1) * K_basis + 1
    end_col <- k * K_basis
    C_integral[k, start_col:end_col] <- c_vec
  }
  
  # Select group-specific coefficients from full coef
  sel_idx <- as.vector(sapply(1:p_total, function(k) (k) + (0:(K_basis - 1)) * X0_dim))
  alpha_hat <- coef[sel_idx]
  V_alpha <- vcov[sel_idx, sel_idx]
  
  # Compute Delta and Cov(Delta)
  delta_hat <- C_integral %*% alpha_hat
  V_delta <- C_integral %*% V_alpha %*% t(C_integral)
  
  # Handle contrast input
  if (is.null(contrast)) {
    C_contrast <- diag(p_total)
  } else {
    C_contrast <- if (is.vector(contrast)) matrix(contrast, nrow = 1) else contrast
    if (ncol(C_contrast) != p_total) stop("Contrast dimension mismatch.")
    if (any(abs(rowSums(C_contrast)) > 1e-8))
      warning("Contrast rows should sum to 0 for interpretability.")
  }
  
  # Apply linear contrast
  theta_hat <- C_contrast %*% delta_hat
  V_theta <- C_contrast %*% V_delta %*% t(C_contrast)
  
  # Compute Wald statistic
  T2 <- t(theta_hat) %*% solve(V_theta) %*% theta_hat
  df <- nrow(C_contrast)
  p_value <- pchisq(T2, df = df, lower.tail = FALSE)
  
  return(list(
    statistic = as.numeric(T2),
    df = df,
    p_value = p_value,
    contrast = C_contrast,
    estimate = theta_hat,
    cov_estimate = V_theta
  ))
}
