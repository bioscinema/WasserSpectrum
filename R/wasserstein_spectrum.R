#' Wasserstein Spectrum Regression via Generalized Least Squares (GLS)
#'
#' Estimates the quantile-varying exposure effect on alpha diversity using functional regression.
#' This method fits a generalized least squares model with basis expansion over quantile levels,
#' allowing for smooth inference on the coefficient function \eqn{\beta_1(t)} with robust standard errors.
#'
#' @param df A data frame containing the alpha diversity values, exposure variable, and optional confounders.
#' @param feature_col A string specifying the column name in \code{df} corresponding to the feature of interest (e.g., taxon, gene, or diversity index).
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
#'                                    feature_col = "Shannon",
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
                                     feature_col, 
                                     outcome_col, 
                                     confounder_cols = NULL,
                                     basis_df = 6, # use < 6 if sample size < 100 to avoid overfitting
                                     t_grid = seq(0.01, 0.99, length.out = 100),
                                     alpha = 0.05,
                                     plot = TRUE,
                                     seed = 123) {
  set.seed(seed)
  
  # Extract data
  y <- df[[feature_col]]
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
      df_lolli <- data.frame(
        quantile = t_grid,
        beta      = beta1_grid,
        lower     = lower,
        upper     = upper
      )
      df_lolli$significant <- with(df_lolli, lower > 0 | upper < 0)
      df_lolli$bar_color <- ifelse(df_lolli$significant, "#e41a1c", "grey70")

      lollipop_p <- ggplot(df_lolli, aes(x = quantile, y = beta)) +
        geom_segment(
          aes(xend = quantile, yend = 0, color = bar_color),
          size = 0.7, lineend = "round", show.legend = FALSE
        ) +
        geom_point(
          shape = 21, fill = "white", color = "#397d54",
          size = 1.5, stroke = 1
        ) +
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          width = 0.01, color = "#397d54", alpha = 0.7
        ) +
        scale_color_identity() +
        labs(
          x = "Quantile level t",
          y = expression(hat(beta)[1](t)),
          title = "Lollipop Plot of Quantile-wise Effect Sizes"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 12,
            margin = margin(b = 8)
          ),
          axis.title.x = element_text(size = 11, margin = margin(t = 6)),
          axis.title.y = element_text(size = 11, margin = margin(r = 6)),
          axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank()
        )
      df_bar <- data.frame(
        angle = t_grid * 360,
        beta  = beta1_grid
      )
      
      # Circular bar plot with custom colors
      circular_p <- ggplot(df_bar, aes(x = angle, y = beta, fill = beta)) +
        geom_col(width = 2, color = NA) +
        coord_polar(theta = "x", start = 0, direction = 1, clip = "off") +
        
        # X axis: quantile as percent
        scale_x_continuous(
          limits = c(0, 360),
          breaks = seq(0, 360, by = 90),
          labels = c("0%", "25%", "50%", "75%", "100%")
        ) +
        
        # Y axis
        scale_y_continuous(
          expand = c(0, 0)
        ) +
        
        # Custom diverging color gradient
        scale_fill_gradient2(
          low = "#d7312d", mid = "#fee395", high = "#6090c1",
          midpoint = 0, name = expression(hat(beta)[1](t))
        ) +
        
        # Labels and theme
        labs(
          title = "Quantile-wise Effect Sizes Across the Distribution",
          x = NULL, y = NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 11, vjust = -1.8),
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 9),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width  = unit(0.3, "cm"),
          plot.margin = margin(10, 20, 10, 20)
        )
    
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
    basis_df = basis_df,
    lollipop_plot = lollipop_p,
    circular_plot = circular_p
  ))
}
