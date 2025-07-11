#' Shape Contrast Test for Wasserstein Spectrum (Multiclass GLS)
#'
#' Performs a Wald-type test for shape-based contrasts of group-specific spectrum curves
#' using a linear combination of quantile levels and group contrasts.
#'
#' @param spectrum_obj An object returned by \code{\link{wasserstein_spectrum_multiclass}},
#' containing basis expansion, estimated coefficients, and covariance matrix.
#' @param t_vec A numeric vector of quantile levels \(t_1, \dots, t_m\) over which to define the shape contrast.
#' @param contrast_t A numeric vector of weights (length \(m\)) defining the shape contrast over quantiles.
#' This vector should typically sum to zero for pure shape comparisons.
#' @param contrast_group Optional group contrast matrix or vector (default is global test across all groups).
#' Each row defines a linear contrast across groups. Rows should sum to 0 for pairwise interpretation.
#'
#' @return A list containing:
#' \item{statistic}{Wald test statistic based on the specified contrast.}
#' \item{df}{Degrees of freedom for the chi-squared test (equal to number of contrast rows).}
#' \item{p_value}{P-value from the chi-squared distribution.}
#' \item{delta_shape}{Contrast-weighted effect estimates across all groups.}
#' \item{cov_shape}{Covariance matrix of \code{delta_shape}.}
#' \item{contrast_group}{The contrast matrix used to compare groups.}
#' \item{estimate}{Final linear contrast estimate (e.g., group difference in shape).}
#' \item{cov_estimate}{Covariance matrix of the final contrast estimate.}
#'
#' @details
#' This function enables formal comparison of the *shape* of spectrum curves (e.g., U-shape vs. linear shift)
#' across groups using user-defined weights over quantiles (\code{contrast_t}) and across groups
#' (\code{contrast_group}). The test statistic is computed as a Wald quadratic form:
#' \deqn{
#' \theta^\top \mathrm{Var}^{-1}(\theta) \theta,
#' }
#' where \(\theta = C_g C_t^\top \Phi(t)\hat\beta\).
#'
#' If \code{contrast_group} is not provided, a global test across all groups is conducted.
#'
#' @examples
#' \dontrun{
#' # U-shaped contrast over quantiles
#' t_vec <- seq(0.1, 0.9, length.out = 5)
#' contrast_t <- c(1, -0.5, -1, -0.5, 1)
#'
#' # Test difference in U-shape between two groups
#' contrast_group <- c(1, -1, 0)  # group 1 vs group 2
#'
#' result <- shape_contrast_test(spectrum_obj = fit,
#'                               t_vec = t_vec,
#'                               contrast_t = contrast_t,
#'                               contrast_group = contrast_group)
#' print(result$p_value)
#' }
#'
#' @importFrom stats pchisq
#' @export
shape_contrast_test <- function(spectrum_obj, t_vec, contrast_t, contrast_group = NULL) {
  Phi <- spectrum_obj$basis
  coef <- spectrum_obj$coef
  vcov <- spectrum_obj$vcov
  comparison_levels <- spectrum_obj$comparison_levels
  
  K_group <- length(comparison_levels) + 1
  p_total <- K_group - 1
  K_basis <- ncol(Phi)
  X0_dim <- length(coef) / K_basis
  
  # Construct shape contrast vector c_final = C_t^T %*% Phi(t_vec)
  if (abs(sum(contrast_t)) > 1e-8) warning("contrast_t should sum to 0 for shape-based tests")
  
  Phi_t <- apply(matrix(t_vec), 1, function(t) as.numeric(predict(Phi, newx = t)))  # K_basis x m
  c_final <- as.vector(contrast_t %*% t(Phi_t))  # shape contrast: 1 x K_basis
  
  # Build shape contrast block matrix
  C_shape <- matrix(0, nrow = p_total, ncol = p_total * K_basis)
  for (k in 1:p_total) {
    start <- (k - 1) * K_basis + 1
    end <- k * K_basis
    C_shape[k, start:end] <- c_final
  }
  
  # Extract alpha and V_alpha
  sel_idx <- as.vector(sapply(1:p_total, function(k) (k) + (0:(K_basis - 1)) * X0_dim))
  alpha_hat <- coef[sel_idx]
  V_alpha <- vcov[sel_idx, sel_idx]
  
  delta_shape <- C_shape %*% alpha_hat
  V_shape <- C_shape %*% V_alpha %*% t(C_shape)
  
  # Global vs contrast test
  if (is.null(contrast_group)) {
    Cg <- diag(p_total)  # Global test: all shape deltas == 0
  } else {
    Cg <- if (is.vector(contrast_group)) matrix(contrast_group, nrow = 1) else contrast_group
    if (ncol(Cg) != p_total) stop("contrast_group has incorrect dimension")
    if (any(abs(rowSums(Cg)) > 1e-8)) warning("contrast_group rows should sum to 0 for interpretability")
  }
  
  theta_hat <- Cg %*% delta_shape
  V_theta <- Cg %*% V_shape %*% t(Cg)
  
  stat <- as.numeric(t(theta_hat) %*% solve(V_theta) %*% theta_hat)
  df <- nrow(Cg)
  p_val <- pchisq(stat, df = df, lower.tail = FALSE)
  
  return(list(
    statistic = stat,
    df = df,
    p_value = p_val,
    delta_shape = delta_shape,
    cov_shape = V_shape,
    contrast_group = Cg,
    estimate = theta_hat,
    cov_estimate = V_theta
  ))
}
