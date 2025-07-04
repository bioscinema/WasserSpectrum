#' Wasserstein Spectrum Regression via Generalized Least Squares (GLS)
#'
#' Estimates the quantile-varying exposure effect on alpha diversity using functional regression.
#' This method fits a generalized least squares model with basis expansion over quantile levels,
#' allowing for smooth inference on the coefficient function \eqn{\beta_1(t)} with robust standard errors.
#'
#' @param df A data frame containing the alpha diversity values, exposure variable, and optional confounders.
#' @param diversity_col A string specifying the column name containing the continuous alpha diversity index.
#' @param outcome_col A string specifying the name of the exposure variable (binary or numeric).
#' @param confounder_cols Optional character vector of column names to adjust for as covariates (default is \code{NULL}).
#' @param basis_df Degrees of freedom used for the B-spline basis in quantile expansion (default is 6).
#' Use smaller values (e.g., 4–5) when sample size is small to avoid overfitting.
#' @param t_grid Numeric vector of quantile levels \eqn{t \in (0,1)} to evaluate along the diversity distribution
#' (default is 100 equally spaced points between 0.01 and 0.99).
#' @param alpha Significance level for confidence intervals (default is 0.05).
#' @param plot Logical. If \code{TRUE}, plots the estimated coefficient curve and confidence bands (default is \code{TRUE}).
#' @param seed Random seed for reproducibility (default is 123).
#'
#' @return A list with the following components:
#' \item{quantiles}{Vector of quantile levels \(t\) used in the regression.}
#' \item{beta1}{Estimated exposure effect \(\hat{\beta}_1(t)\) across the quantile spectrum.}
#' \item{se1}{Robust standard errors for \(\hat{\beta}_1(t)\).}
#' \item{lower}{Lower bound of the \((1 - \alpha)\)-level confidence interval.}
#' \item{upper}{Upper bound of the \((1 - \alpha)\)-level confidence interval.}
#' \item{basis}{The B-spline basis matrix used for quantile expansion.}
#' \item{coef}{Full GLS coefficient vector for all basis-covariate combinations.}
#' \item{vcov}{Robust sandwich variance-covariance matrix for the estimated coefficients.}
#' \item{idx_x}{Index of the exposure variable in the model matrix.}
#' \item{basis_df}{Degrees of freedom used for the spline basis.}
#'
#' @details
#' The function approximates the conditional quantile function of diversity using a B-spline basis
#' over quantile levels and fits a linear model to binary indicators of whether diversity exceeds
#' each quantile threshold. The resulting effect curve \(\hat{\beta}_1(t)\) characterizes how an
#' exposure influences different portions of the diversity distribution.
#'
#' Robust standard errors are calculated using heteroskedasticity-consistent estimators to account
#' for residual correlation across quantile levels.
#'
#' @examples
#' \dontrun{
#' result <- wasserstein_spectrum(df = example_df,
#'                                    diversity_col = "Shannon",
#'                                    outcome_col = "BMI",
#'                                    confounder_cols = c("Age", "Sex"))
#' plot(result$quantiles, result$beta1, type = "l")
#' }
#'
#' @importFrom stats quantile model.matrix coef lm qnorm
#' @importFrom splines bs
#' @importFrom sandwich vcovHC
#' @export
wasserstein_spectrum <- function(df, 
                                     diversity_col, 
                                     outcome_col, 
                                     confounder_cols = NULL,
                                     basis_df = 6, # use < 6 if sample size < 100 to avoid overfitting
                                     t_grid = seq(0.01, 0.99, length.out = 100),
                                     alpha = 0.05,
                                     plot = TRUE,
                                     seed = 123) {
  set.seed(seed)
  
  library(splines)
  
  # Extract data
  y <- df[[diversity_col]]
  x <- df[[outcome_col]]
  if (is.factor(x)) {
    if (nlevels(x) != 2) stop("outcome_col must be binary or numeric.")
    x <- as.numeric(x) - 1
  } else if (!is.numeric(x)) {
    stop("outcome_col must be numeric or a binary factor.")
  }
  conf_df <- if (!is.null(confounder_cols)) df[, confounder_cols, drop = FALSE] else NULL
  
  n <- nrow(df)
  m <- length(t_grid)
  
  # Build quantile thresholds
  quantile_thresholds <- quantile(y, probs = t_grid, type = 8)
  
  # Construct binary matrix: y_bin[i,j] = 1 if y_i > q_j
  y_bin_mat <- outer(y, quantile_thresholds, FUN = ">") * 1  # n x m
  
  # Construct design matrix: for each quantile level, same X
  model_data <- data.frame(x = x)
  if (!is.null(conf_df)) model_data <- cbind(model_data, conf_df)
  X0 <- model.matrix(~ ., data = model_data)  # n x p
  p <- ncol(X0)
  
  # Construct B-spline basis on t_grid
  Phi <- bs(t_grid, df = basis_df, intercept = TRUE)  # m x K
  K <- ncol(Phi)
  
  # Build full model matrix: kronecker product for beta1(t)
  X_big <- kronecker(Phi, X0)  # (n*m) x (p*K)
  y_vec <- as.vector(y_bin_mat)  # flatten row-wise (n*m)
  
  # Fit model
  fit <- lm(y_vec ~ X_big - 1)
  coef_hat <- coef(fit)  # length p*K
  vcov_hat <- sandwich::vcovHC(fit, type = "HC1")
  
  # Extract beta1(t): second column of X0 × basis
  idx_x <- which(colnames(X0) == "x")
  if (length(idx_x) == 0) stop("Couldn't find predictor `x` in design matrix.")
  
  beta1_grid <- numeric(m)
  se1_grid <- numeric(m)
  
  for (j in 1:m) {
    phi_j <- Phi[j, ]  # 1 x K
    # Construct selector vector: 1 at positions (idx_x + (0:(K-1)) * p)
    sel_idx <- idx_x + (0:(K-1)) * p
    beta1_grid[j] <- sum(coef_hat[sel_idx] * phi_j)
    se1_grid[j] <- sqrt(t(phi_j) %*% vcov_hat[sel_idx, sel_idx] %*% phi_j)
  }
  
  lower <- beta1_grid - qnorm(1 - alpha / 2) * se1_grid
  upper <- beta1_grid + qnorm(1 - alpha / 2) * se1_grid
  
  # Plot
  if (plot) {
    plot(t_grid, beta1_grid, type = "l", lwd = 2, col = "darkred",
         ylim = range(c(lower, upper)), xlab = "Quantile level t",
         ylab = expression(beta[1](t)), main = "Wasserstein Spectrum (GLS)")
    polygon(c(t_grid, rev(t_grid)),
            c(upper, rev(lower)),
            col = adjustcolor("darkred", alpha.f = 0.2), border = NA)
    abline(h = 0, lty = 2)
    legend("topright", legend = c(expression(hat(beta)[1](t)), "CI"),
           col = c("darkred", adjustcolor("darkred", alpha.f = 0.2)), lwd = c(2, 8), bty = "n")
  }
  
  return(list(
    quantiles = t_grid,
    beta1 = beta1_grid,
    se1 = se1_grid,
    lower = lower,
    upper = upper,
    basis = Phi,
    coef = coef_hat,
    vcov = vcov_hat,
    idx_x = idx_x,
    basis_df = basis_df
  ))
}
