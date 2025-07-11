#' Interval-wise MANOVA Test for Group Differences in Spectrum Curves
#'
#' Performs a Wald-type global MANOVA test to assess whether the integrated spectrum curves differ
#' across groups (vs reference) over a specified quantile interval \eqn{[a, b]}.
#'
#' @param spectrum_obj An object returned by \code{\link{wasserstein_spectrum_multiclass}}.
#' It must contain estimated coefficients, basis functions, quantile grid, and variance-covariance matrix.
#' @param a Lower bound of the quantile interval \eqn{a \in (0, 1)}.
#' @param b Upper bound of the quantile interval \eqn{b > a}.
#'
#' @return A list with the following components:
#' \item{statistic}{Wald statistic (Hotelling-type \eqn{T^2}) for the global test.}
#' \item{df}{Degrees of freedom used in the test.}
#' \item{p_value}{P-value under the chi-squared distribution.}
#' \item{delta}{Integrated effect estimates (vector of group vs reference differences over \eqn{[a, b]}).}
#' \item{cov_delta}{Covariance matrix of \code{delta}.}
#'
#' @details
#' This function tests whether group-specific spectrum curves differ significantly from the reference group
#' when integrated over the interval \eqn{[a, b]}. It uses spline basis projection and multivariate Wald-type
#' inference.
#'
#' The null hypothesis is:
#' \deqn{H_0: \int_a^b \beta_k(t)\,dt = 0 \quad \text{for all } k,}
#' where \eqn{\beta_k(t)} denotes the deviation curve for group \eqn{k} vs reference.
#'
#' @seealso \code{\link{if_manova_contrast}} for testing user-specified contrasts instead of the global null.
#'
#' @examples
#' \dontrun{
#' fit <- wasserstein_spectrum_multiclass(...)
#' test <- if_manova(fit, a = 0.2, b = 0.8)
#' print(test$p_value)
#' }
#'
#' @importFrom stats pchisq
#' @export
if_manova <- function(spectrum_obj, a, b) {
  # Extract components
  Phi <- spectrum_obj$basis          # m x K (basis matrix)
  t_grid <- spectrum_obj$quantiles   # length m
  coef <- spectrum_obj$coef          # length p*K
  vcov <- spectrum_obj$vcov          # (p*K) x (p*K)
  comparison_levels <- spectrum_obj$comparison_levels
  K_group <- length(comparison_levels) + 1  # including reference group
  K_basis <- ncol(Phi)              # number of basis functions
  p_total <- K_group - 1            # number of effect curves vs reference
  X0_dim <- length(coef) / K_basis  # original model matrix columns
  
  # Build contrast matrix C: (K-1) x ((K-1)*K_basis)
  idx_ab <- which(t_grid >= a & t_grid <= b)
  Phi_sub <- Phi[idx_ab, , drop = FALSE]  # subset of basis
  delta_t <- mean(diff(t_grid))
  c_vec <- colSums(Phi_sub) * delta_t  # numeric integration
  
  # Create block matrix: each row picks coefficients of one group
  C_mat <- matrix(0, nrow = p_total, ncol = p_total * K_basis)
  for (k in 1:p_total) {
    start_col <- (k - 1) * K_basis + 1
    end_col <- k * K_basis
    C_mat[k, start_col:end_col] <- c_vec
  }
  
  # Select alpha_hat and Cov_alpha from full coef and vcov
  sel_idx <- as.vector(sapply(1:p_total, function(k) {
    (k) + (0:(K_basis - 1)) * X0_dim
  }))
  alpha_hat <- coef[sel_idx]
  V_alpha <- vcov[sel_idx, sel_idx]
  
  # Compute delta and its covariance
  delta_hat <- as.vector(C_mat %*% alpha_hat)
  V_delta <- C_mat %*% V_alpha %*% t(C_mat)
  
  # Wald-type statistic
  T2 <- t(delta_hat) %*% solve(V_delta) %*% delta_hat
  df <- p_total
  p_value <- pchisq(T2, df = df, lower.tail = FALSE)
  
  return(list(
    statistic = as.numeric(T2),
    df = df,
    p_value = p_value,
    delta = delta_hat,
    cov_delta = V_delta
  ))
}
