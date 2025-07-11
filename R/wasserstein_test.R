#' Permutation-Based Wasserstein Test for Two-Group Alpha Diversity Comparison
#'
#' Performs a two-sample hypothesis test for distributional differences in alpha diversity
#' between two groups using the 1D Wasserstein distance (also known as Earth Mover’s Distance).
#' This test is robust to ties and sensitive to shifts in location, scale, and shape of the
#' diversity distribution. The null distribution is generated via permutation of group labels.
#'
#' @param df A data frame containing the alpha diversity values and group labels.
#' @param feature_col A string specifying the column name in \code{df} corresponding to the feature of interest (e.g., taxon, gene, or diversity index).
#' @param group_col A string specifying the column name in \code{df} indicating the group variable.
#' Must have exactly two unique values.
#' @param nperm Number of permutations used to generate the null distribution (default is 1000).
#' @param plot Logical indicating whether to display a histogram of the null distribution
#' with the observed statistic marked (default is \code{TRUE}).
#' @param seed An integer seed for reproducibility (default is 123).
#'
#' @return A list with the following components:
#' \item{statistic}{Observed Wasserstein distance between the two groups.}
#' \item{p_value}{Permutation p-value comparing observed statistic to null distribution.}
#' \item{null_distribution}{A numeric vector containing the permuted test statistics.}
#'
#' @details
#' The test computes the 1D Wasserstein distance between the empirical distributions
#' of the alpha diversity values in the two groups. A permutation-based p-value is calculated
#' by comparing the observed distance to the distribution of distances obtained from permuted group labels.
#'
#' This approach is particularly well-suited to microbiome studies, where diversity distributions
#' may exhibit heavy ties, skewness, and heterogeneity that are poorly captured by traditional
#' location-based tests.
#'
#' @examples
#' \dontrun{
#' result <- wasserstein_test(df = example_df, feature_col = "Shannon", group_col = "Group")
#' print(result$p_value)
#' }
#'
#' @importFrom transport wasserstein1d
#' @export
wasserstein_test <- function(df, feature_col, group_col, 
                             nperm = 1000, plot = TRUE, seed = 123) {
  set.seed(seed)
  
  # Extract variables
  diversity <- df[[feature_col]]
  group <- df[[group_col]]
  groups <- unique(group)
  if (length(groups) != 2) stop("group_col must have exactly 2 groups")
  
  x <- diversity[group == groups[1]]
  y <- diversity[group == groups[2]]
  
  # Observed Wasserstein distance
  T_obs <- transport::wasserstein1d(x, y)
  
  # Permutation null distribution
  all_vals <- c(x, y)
  n_x <- length(x)
  null_dist <- numeric(nperm)
  for (b in seq_len(nperm)) {
    perm <- sample(all_vals)
    perm_x <- perm[1:n_x]
    perm_y <- perm[(n_x + 1):length(perm)]
    null_dist[b] <- transport::wasserstein1d(perm_x, perm_y)
  }
  
  # P-value
  p_val <- mean(null_dist >= T_obs)
  
  # Plot
  if (plot) {
    hist(null_dist, breaks = 40, col = "lightblue", border = "white",
         main = "Null distribution of Wasserstein statistic",
         xlab = "Wasserstein distance")
    abline(v = T_obs, col = "red", lwd = 2)
    legend("topright", legend = c("Observed"), col = "red", lwd = 2, bty = "n")
  }
  
  return(list(
    statistic = T_obs,
    p_value = p_val,
    null_distribution = null_dist
  ))
}