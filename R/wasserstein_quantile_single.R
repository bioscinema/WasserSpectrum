#' Wasserstein Quantile-Based Inference for a Single Feature
#'
#' This function performs quantile-based inference on a single feature (e.g., taxon, gene, or metabolite)
#' using the Wasserstein spectrum framework. It estimates the association between an outcome variable and
#' a feature’s abundance across quantiles via B-spline projection, and aggregates quantile-level p-values
#' using the (adaptive) Cauchy combination method. It also returns the area under the estimated effect curve.
#'
#' @param df A data frame or matrix with features as columns and samples as rows.
#' @param feature_col Character. The column name in `df` corresponding to the feature of interest.
#' @param metadata A data frame containing metadata with sample-level variables.
#' @param outcome_col Character. The column name in `metadata` for the primary outcome.
#' @param confounder_cols Optional character vector. Column names in `metadata` for confounder adjustment.
#' @param basis_df Integer. Degrees of freedom for the B-spline basis (default = 6).
#' @param t_grid Numeric vector. Quantile levels to evaluate (default = seq(0.01, 0.99, length.out = 100)).
#' @param adaptive Logical. Whether to use adaptive Cauchy combination (default = FALSE).
#' @return A data frame with columns: feature name, global p-value, AUC, and absolute AUC.
#' @importFrom splines bs
#' @importFrom sandwich vcovHC
#' @importFrom pracma trapz
#' @export
wasserstein_quantile_single <- function(df,
                                        feature_col,
                                        metadata,
                                        outcome_col,
                                        confounder_cols = NULL,
                                        basis_df = 6,
                                        t_grid = seq(0.01, 0.99, length.out = 100),
                                        adaptive = FALSE) {
  requireNamespace("splines")
  requireNamespace("sandwich")
  requireNamespace("pracma")
  
  Phi <- splines::bs(t_grid, df = basis_df, intercept = TRUE)
  K <- ncol(Phi)
  m <- length(t_grid)
  
  y <- df[[feature_col]]
  x <- metadata[[outcome_col]]
  
  if (is.factor(x)) {
    if (nlevels(x) != 2) stop("Factor outcome must have 2 levels.")
    x <- as.numeric(x) - 1
  } else if (!is.numeric(x)) {
    stop("Outcome must be numeric or binary factor.")
  }
  
  conf_df <- if (!is.null(confounder_cols)) metadata[, confounder_cols, drop = FALSE] else NULL
  
  quantile_thresholds <- quantile(y, probs = t_grid, type = 8)
  y_bin_mat <- outer(y, quantile_thresholds, FUN = ">") * 1
  
  model_data <- data.frame(x = x)
  if (!is.null(conf_df)) model_data <- cbind(model_data, conf_df)
  X0 <- model.matrix(~ ., data = model_data)
  p_dim <- ncol(X0)
  
  X_big <- kronecker(Phi, X0)
  y_vec <- as.vector(y_bin_mat)
  
  fit <- lm(y_vec ~ X_big - 1)
  coef_hat <- coef(fit)
  vcov_hat <- sandwich::vcovHC(fit, type = "HC1")
  
  idx_x <- which(colnames(X0) == "x")
  if (length(idx_x) == 0) stop("Outcome variable 'x' not found.")
  sel_idx <- idx_x + (0:(K - 1)) * p_dim
  
  beta1 <- as.numeric(Phi %*% coef_hat[sel_idx])
  se1 <- sapply(1:m, function(j) {
    phi_j <- Phi[j, , drop = FALSE]
    sqrt(phi_j %*% vcov_hat[sel_idx, sel_idx] %*% t(phi_j))
  })
  
  z_scores <- beta1 / se1
  pvals <- 2 * (1 - pnorm(abs(z_scores)))
  
  if (!adaptive) {
    weights <- rep(1 / m, m)
    cauchy_stat <- sum(weights * tan((0.5 - pvals) * pi))
    p_cauchy <- 0.5 - atan(cauchy_stat) / pi
  } else {
    k_list <- floor(m * c(0.05, 0.1, 0.2, 0.3, 0.5))
    p_subset <- sapply(k_list, function(k) {
      top_k <- sort(pvals)[1:k]
      stat_k <- sum(tan((0.5 - top_k) * pi)) / k
      0.5 - atan(stat_k) / pi
    })
    stat2 <- sum(tan((0.5 - p_subset) * pi)) / length(p_subset)
    p_cauchy <- 0.5 - atan(stat2) / pi
  }
  
  auc <- pracma::trapz(t_grid, beta1)
  abs_auc <- pracma::trapz(t_grid, abs(beta1))
  
  data.frame(
    Feature = feature_col,
    p_value = p_cauchy,
    auc = auc,
    abs_auc = abs_auc
  )
}
