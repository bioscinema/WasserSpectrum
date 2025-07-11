#' Spectrum-Based MANOVA (Hotelling-Type Test) for Multiclass GLS Regression
#'
#' Performs a Hotelling-type test at each quantile to assess whether group-level contrast curves
#' \eqn{\beta_k(t)} jointly differ from zero. This function applies to results from
#' \code{\link{wasserstein_spectrum_multiclass}}.
#'
#' @param spectrum_obj An object returned by \code{\link{wasserstein_spectrum_multiclass}},
#' containing estimated contrast curves, robust covariance matrix, and spline basis matrix.
#'
#' @return A data frame with the following columns:
#' \item{quantile}{The quantile level \(t\) at which the test is performed.}
#' \item{Fval}{The Hotelling-type test statistic at each quantile.}
#' \item{pval}{The corresponding p-value under an F-distribution with \code{df1 = K - 1} and \code{df2 = Inf}.}
#'
#' @details
#' For each quantile level \eqn{t_j}, the function tests the null hypothesis that
#' all group-level contrast effects \eqn{\beta_k(t_j)} are simultaneously zero.
#' The test statistic is a scaled version of Hotelling’s \(T^2\), adapted for robust covariance:
#' \deqn{
#' F(t_j) = \frac{1}{K - 1} \beta(t_j)^\top \Sigma^{-1}(t_j) \beta(t_j),
#' }
#' where \(\Sigma(t_j)\) is the estimated covariance matrix of \(\beta(t_j)\).
#'
#' This procedure yields a quantile-wise map of where multigroup differences are strongest.
#'
#' @examples
#' \dontrun{
#' fit <- wasserstein_spectrum_multiclass(
#'   df = example_df,
#'   feature_col = "Shannon",
#'   outcome_col = "Group",
#'   reference_level = "Healthy",
#'   confounder_cols = c("Age", "Sex")
#' )
#' manova_result <- spectrum_manova(fit)
#' plot(manova_result$quantile, -log10(manova_result$pval), type = "l")
#' }
#'
#' @importFrom stats pf
#' @export
spectrum_manova <- function(spectrum_obj) {
  beta <- spectrum_obj$beta         # (K-1) x m
  vcov <- spectrum_obj$vcov         # full (p*K) x (p*K)
  coef <- spectrum_obj$coef         # vectorized length p*K
  basis <- spectrum_obj$basis       # m x K_basis
  t_grid <- spectrum_obj$quantiles
  comparison_levels <- spectrum_obj$comparison_levels
  
  m <- length(t_grid)
  K_basis <- ncol(basis)
  p_total <- length(comparison_levels)       # number of groups vs reference
  n_coef <- length(coef)
  p <- n_coef / K_basis                      # total number of columns in X0
  
  # Which positions in coef correspond to group-specific predictors?
  # Based on your code, they are always (k + 1) + (0:(K_basis - 1)) * p, for k = 1 to p_total
  group_sel_mat <- sapply(1:p_total, function(k) {
    (k + 1) + (0:(K_basis - 1)) * p
  })  # a matrix: K_basis x p_total
  
  Fvals <- numeric(m)
  pvals <- numeric(m)
  
  for (j in 1:m) {
    phi_j <- basis[j, ]  # 1 x K_basis
    beta_j <- beta[, j]  # vector of length p_total
    
    Sigma_j <- matrix(0, p_total, p_total)
    for (k in 1:p_total) {
      sel_k <- group_sel_mat[, k]
      Sigma_kj <- as.numeric(t(phi_j) %*% vcov[sel_k, sel_k] %*% phi_j)
      Sigma_j[k, k] <- Sigma_kj
    }
    
    # F-type test (Hotelling's T^2)
    Fval <- t(beta_j) %*% solve(Sigma_j) %*% beta_j / p_total
    pval <- 1 - pf(Fval, df1 = p_total, df2 = Inf)
    Fvals[j] <- Fval
    pvals[j] <- pval
  }
  
  return(data.frame(
    quantile = t_grid,
    Fval = Fvals,
    pval = pvals
  ))
}
