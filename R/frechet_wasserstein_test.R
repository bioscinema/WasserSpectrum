#' Fréchet–Wasserstein Test for Multi-Group Alpha Diversity Comparison
#'
#' Performs a permutation-based test for comparing alpha diversity distributions across
#' multiple groups using the Fréchet variance in Wasserstein space. This method generalizes
#' the two-sample Wasserstein test to \(G > 2\) groups by quantifying the group-wise
#' dispersion around a shared Wasserstein barycenter.
#'
#' @param df A data frame containing the alpha diversity values and group labels.
#' @param diversity_col A string specifying the column name in \code{df} that contains the alpha diversity index.
#' @param group_col A string specifying the column name in \code{df} that indicates the group variable. Must have at least two unique levels.
#' @param nperm Number of permutations used to generate the null distribution (default is 1000).
#' @param plot Logical. If \code{TRUE}, displays a histogram of the null distribution with the observed test statistic (default is \code{TRUE}).
#' @param seed An integer to set the random seed for reproducibility (default is 123).
#'
#' @return A list with the following components:
#' \item{T_obs}{Observed Fréchet–Wasserstein test statistic, measuring total group-wise dispersion.}
#' \item{p_value}{Permutation p-value comparing \code{T_obs} to the null distribution.}
#' \item{T_null}{A numeric vector containing the permuted test statistics under the null hypothesis.}
#'
#' @details
#' The Fréchet–Wasserstein test extends classical ANOVA to the space of distributions by computing
#' the squared Wasserstein distances between group-wise quantile functions and their shared barycenter.
#' The observed statistic is compared to a null distribution obtained by permuting group labels,
#' providing a robust and distribution-free inference procedure. This method is particularly suited
#' to microbiome studies, where diversity values often exhibit ties, skewness, or non-Gaussian behavior.
#'
#' @examples
#' \dontrun{
#' result <- frechet_wasserstein_test(df = example_df,
#'                                    diversity_col = "Shannon",
#'                                    group_col = "Group")
#' print(result$p_value)
#' }
#'
#' @importFrom stats quantile
#' @export
frechet_wasserstein_test <- function(df, 
                                     diversity_col, 
                                     group_col, 
                                     nperm = 1000, 
                                     plot = TRUE,
                                     seed = 123) {
  set.seed(seed)
  library(dplyr)
  
  y <- df[[diversity_col]]
  group <- as.factor(df[[group_col]])
  groups <- levels(group)
  G <- length(groups)
  n <- length(y)
  
  # quantile grid
  t_grid <- seq(0.01, 0.99, length.out = 100)
  
  # quantile functions
  group_quantiles <- lapply(groups, function(g) {
    y_g <- y[group == g]
    quantile(y_g, probs = t_grid, type = 8)
  })
  names(group_quantiles) <- groups
  
  # Fréchet barycenter
  q_bar <- Reduce("+", group_quantiles) / G
  
  # Wasserstein² distance for each group
  D_gs <- sapply(group_quantiles, function(q_g) mean((q_g - q_bar)^2))
  n_gs <- table(group)
  T_obs <- sum(n_gs * D_gs)
  
  # Permutation null distribution
  T_null <- numeric(nperm)
  for (b in 1:nperm) {
    group_perm <- sample(group)
    perm_qs <- lapply(groups, function(g) {
      y_g <- y[group_perm == g]
      quantile(y_g, probs = t_grid, type = 8)
    })
    q_bar_perm <- Reduce("+", perm_qs) / G
    D_gs_perm <- sapply(perm_qs, function(q_g) mean((q_g - q_bar_perm)^2))
    T_null[b] <- sum(n_gs * D_gs_perm)
  }
  
  # P-value
  p_val <- mean(T_null >= T_obs)
  
  # Plot
  if (plot) {
    hist(T_null, breaks = 40, col = "gray80", border = "white",
         main = "Fréchet Wasserstein Test",
         xlab = "Test Statistic under Null")
    abline(v = T_obs, col = "red", lwd = 2)
    legend("topright", legend = paste0("Observed = ", round(T_obs, 4),
                                       "\nP = ", signif(p_val, 4)),
           bty = "n")
  }
  
  return(list(
    T_obs = T_obs,
    p_value = p_val,
    T_null = T_null
  ))
}
