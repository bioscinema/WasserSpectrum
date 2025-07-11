#' Multiclass Wasserstein Spectrum Regression via Generalized Least Squares (GLS)
#'
#' Estimates the quantile-resolved effect of a categorical exposure variable on alpha diversity
#' using functional regression with basis expansion and GLS. The model captures group-specific
#' contrast curves \eqn{\beta_k(t)} comparing each level to the user-defined reference group.
#'
#' @param df A data frame containing the diversity index, exposure variable, and optional confounders.
#' @param feature_col A string specifying the column name in \code{df} corresponding to the feature of interest (e.g., taxon, gene, or diversity index).
#' @param outcome_col A string specifying the name of the categorical exposure variable (with \code{K >= 2} levels).
#' @param reference_level A string indicating which level of the exposure variable to use as the reference group.
#' @param confounder_cols Optional character vector of column names to include as covariates (default is \code{NULL}).
#' @param basis_df Degrees of freedom for the B-spline basis expansion across quantiles (default is 6).
#' @param t_grid Numeric vector of quantile levels \eqn{t \in (0, 1)} to evaluate (default is 100 points from 0.01 to 0.99).
#' @param alpha Significance level for constructing confidence intervals (default is 0.05).
#' @param plot Logical. If \code{TRUE}, returns a faceted ggplot visualizing the contrast curves (default is \code{TRUE}).
#' @param seed Integer seed for reproducibility (default is 123).
#'
#' @return A list containing:
#' \item{quantiles}{Quantile levels used in the regression.}
#' \item{beta}{A matrix of estimated contrast curves \(\hat{\beta}_k(t)\) for each comparison group (rows) across quantiles (columns).}
#' \item{se}{Matrix of standard errors for the estimated curves.}
#' \item{lower}{Lower bounds of the \((1 - \alpha)\)-level confidence intervals.}
#' \item{upper}{Upper bounds of the \((1 - \alpha)\)-level confidence intervals.}
#' \item{comparison_levels}{Names of the non-reference levels being contrasted.}
#' \item{reference_level}{The reference level used in contrast coding.}
#' \item{basis}{B-spline basis matrix used for quantile expansion.}
#' \item{coef}{Estimated coefficients from the GLS fit.}
#' \item{vcov}{Robust variance-covariance matrix of the estimated coefficients.}
#' \item{plot}{A \code{ggplot} object showing the spectrum plot if \code{plot = TRUE}, otherwise \code{NULL}.}
#'
#' @details
#' This function generalizes the binary GLS spectrum regression to multiclass exposures. At each quantile threshold,
#' binary outcomes are defined based on whether the diversity index falls below the threshold.
#' B-spline basis expansion is used to parameterize smooth group-level effects \(\beta_k(t)\), and
#' coefficient curves are estimated using generalized least squares with robust (sandwich) standard errors.
#'
#' The output enables inference and visualization of how different groups differ across the spectrum of diversity.
#'
#' @examples
#' \dontrun{
#' result <- wasserstein_spectrum_multiclass(
#'   df = example_df,
#'   feature_col = "Shannon",
#'   outcome_col = "DiseaseGroup",
#'   reference_level = "Healthy",
#'   confounder_cols = c("Age", "Sex")
#' )
#' result$plot
#' }
#'
#' @importFrom stats quantile model.matrix coef lm qnorm relevel
#' @importFrom splines bs
#' @importFrom sandwich vcovHC
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline facet_wrap labs theme_minimal
#' @export
wasserstein_spectrum_multiclass <- function(df,
                                                feature_col,
                                                outcome_col,
                                                reference_level,
                                                confounder_cols = NULL,
                                                basis_df = 6,
                                                t_grid = seq(0.01, 0.99, length.out = 100),
                                                alpha = 0.05,
                                                plot = TRUE,
                                                seed = 123) {
  set.seed(seed)
  
  # Outcome and predictors
  y <- df[[feature_col]]
  x <- factor(df[[outcome_col]])
  x <- relevel(x, ref = reference_level)
  x_df <- data.frame(x = x)
  X0 <- model.matrix(~ x, data = x_df)  # includes intercept + K−1 dummies
  coef_names <- colnames(X0)[-1]
  
  if (!is.null(confounder_cols)) {
    conf_df <- df[, confounder_cols, drop = FALSE]
    X0 <- cbind(X0, model.matrix(~ ., data = conf_df)[, -1])
  }
  
  p <- ncol(X0)
  n <- nrow(df)
  m <- length(t_grid)
  
  # Binary outcomes across quantiles
  quantile_thresholds <- quantile(y, probs = t_grid, type = 8)
  y_bin_mat <- outer(y, quantile_thresholds, FUN = ">") * 1  # n x m
  
  # Spline basis on quantile grid
  Phi <- bs(t_grid, df = basis_df, intercept = TRUE)  # m x K
  K <- ncol(Phi)
  
  # Kronecker design matrix: (n*m) x (p*K)
  X_big <- kronecker(Phi, X0)  # (n*m) x (p*K)
  y_vec <- as.vector(y_bin_mat)
  
  # Fit model
  fit <- lm(y_vec ~ X_big - 1)
  coef_hat <- coef(fit)
  vcov_hat <- sandwich::vcovHC(fit, type = "HC1")
  
  # Extract beta_k(t): for each group k (vs reference)
  group_idx <- which(grepl("^x", colnames(model.matrix(~ x))))
  n_comp <- length(group_idx)
  
  beta_mat <- matrix(NA, nrow = n_comp, ncol = m)
  se_mat <- matrix(NA, nrow = n_comp, ncol = m)
  lower_mat <- matrix(NA, nrow = n_comp, ncol = m)
  upper_mat <- matrix(NA, nrow = n_comp, ncol = m)
  
  for (j in 1:m) {
    phi_j <- Phi[j, ]
    for (k in 1:n_comp) {
      idx_k <- group_idx[k]
      sel <- idx_k + (0:(K-1)) * p
      beta_hat <- sum(coef_hat[sel] * phi_j)
      se_hat <- sqrt(t(phi_j) %*% vcov_hat[sel, sel] %*% phi_j)
      beta_mat[k, j] <- beta_hat
      se_mat[k, j] <- se_hat
      lower_mat[k, j] <- beta_hat - qnorm(1 - alpha / 2) * se_hat
      upper_mat[k, j] <- beta_hat + qnorm(1 - alpha / 2) * se_hat
    }
  }
  
  # Plot
  if (plot) {
    plot_df <- data.frame()
    for (k in 1:n_comp) {
      temp <- data.frame(
        quantile = t_grid,
        beta = beta_mat[k, ],
        lower = lower_mat[k, ],
        upper = upper_mat[k, ],
        group = rep(coef_names[k], m)
      )
      plot_df <- rbind(plot_df, temp)
    }
    
    plot_df$group <- factor(plot_df$group, levels = coef_names)
    
    p <- ggplot(plot_df, aes(x = quantile, y = beta)) +
      geom_line(color = "blue", linewidth = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ group, ncol = 1, scales = "free_y") +
      labs(title = "Multiclass Wasserstein Spectrum (GLS)",
           x = "Quantile level t", y = expression(beta[1](t))) +
      theme_minimal(base_size = 14)
  } else {
    p <- NULL
  }
  
  return(list(
    quantiles = t_grid,
    beta = beta_mat,
    se = se_mat,
    lower = lower_mat,
    upper = upper_mat,
    comparison_levels = coef_names,
    reference_level = reference_level,
    basis = Phi,
    coef = coef_hat,
    vcov = vcov_hat,
    plot = p
  ))
}
