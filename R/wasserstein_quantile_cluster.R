#' @useDynLib WasserSpectrum, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom splines bs
#' @importFrom sandwich vcovHC
#' @importFrom stats pnorm quantile model.matrix lm coef residuals p.adjust
#' @importFrom cluster pam
#' @importFrom fpc pamk
#' @importFrom pbapply pblapply
NULL

#' Wasserstein Quantile-Based Inference and Clustering for Microbiome Data
#'
#' This function performs quantile-based inference on each taxon using the Wasserstein spectrum framework.
#' It estimates the effect of a continuous or binary outcome on taxon abundance across quantiles using basis projection,
#' computes Wald-type test statistics, and aggregates p-values via (adaptive) Cauchy combination.
#' Optionally, it clusters significant taxa based on the similarity of their effect curves.
#' The function leverages C++ acceleration for efficiency.
#'
#' @param df A matrix or data frame of taxa (columns) by samples (rows). Transpose if taxa are rows.
#' @param metadata A data frame of sample-level metadata. Must include `outcome_col` and optionally `confounder_cols`.
#' @param outcome_col Character. Column name in `metadata` corresponding to the primary outcome variable.
#' @param confounder_cols Optional character vector of column names in `metadata` for confounder adjustment.
#' @param basis_df Integer. Degrees of freedom for the B-spline basis (default = 6).
#' @param t_grid Numeric vector of quantile levels to evaluate (default = 100 values from 0.01 to 0.99).
#' @param adaptive Logical. If TRUE, applies adaptive Cauchy combination (default = FALSE).
#' @param cluster Logical. If TRUE, clusters significant taxa using L1 distance and PAM (default = FALSE).
#' @param seed Integer for reproducibility (default = 123).
#'
#' @return A data frame with one row per taxon and columns:
#' \describe{
#'   \item{OTU}{The taxon name (i.e., column name from `df`).}
#'   \item{p_value}{Cauchy-combined global p-value testing association with the outcome.}
#'   \item{adj_p_value}{Bonferroni-adjusted p-value.}
#'   \item{Cluster}{Cluster assignment (if `cluster = TRUE`).}
#'   \item{IsMedoid}{TRUE if taxon is a medoid within its cluster (if `cluster = TRUE`).}
#' }
#'
#' @details
#' Internally, for each taxon, a series of binary outcomes is constructed across quantile thresholds,
#' and a linear model is fit with B-spline basis and covariate design via Kronecker product.
#' The coefficient of interest (associated with the outcome) is extracted across quantiles,
#' and test statistics are computed using heteroskedasticity-consistent covariance estimates (HC1).
#'
#' If clustering is enabled, taxa with significant adjusted p-values (< 0.05) are clustered
#' using L1 distance on effect curves and partitioned using `pamk()` from `fpc`.
#'
#' @examples
#' \dontrun{
#' result <- wasserstein_quantile_cluster(
#'   df = t(count_data),
#'   metadata = metadata,
#'   outcome_col = "BMI",
#'   confounder_cols = c("Age", "Gender"),
#'   cluster = TRUE
#' )
#' }
#'
#' @export
wasserstein_quantile_cluster <- function(df,
                                         metadata,
                                         outcome_col,
                                         confounder_cols = NULL,
                                         basis_df = 6,
                                         t_grid = seq(0.01, 0.99, length.out = 100),
                                         adaptive = FALSE,
                                         cluster = FALSE,
                                         seed = 123) {
  set.seed(seed)
  
  otus <- colnames(df)
  Phi <- bs(t_grid, df = basis_df, intercept = TRUE)
  K <- ncol(Phi)
  
  run_for_otu <- function(otu) {
    y <- df[[otu]]
    x <- metadata[[outcome_col]]
    if (is.factor(x)) {
      if (nlevels(x) != 2) return(NULL)
      x <- as.numeric(x) - 1
    } else if (!is.numeric(x)) return(NULL)
    
    conf_df <- if (!is.null(confounder_cols)) metadata[, confounder_cols, drop = FALSE] else NULL
    n <- length(y)
    m <- length(t_grid)
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
    vcov_hat <- vcovHC(fit, type = "HC1")
    
    idx_x <- which(colnames(X0) == "x")
    if (length(idx_x) == 0) return(NULL)
    sel_idx <- idx_x + (0:(K - 1)) * p_dim
    
    res_cpp <- compute_beta1_se(Phi, coef_hat, vcov_hat, sel_idx)
    beta1_grid <- res_cpp$beta1
    se1_grid <- res_cpp$se1
    
    z_scores <- beta1_grid / se1_grid
    pvals <- 2 * (1 - pnorm(abs(z_scores)))
    
    if (!adaptive) {
      weights <- rep(1 / m, m)
      cauchy_stat <- sum(weights * tan((0.5 - pvals) * pi))
      p_cauchy <- 0.5 - atan(cauchy_stat) / pi
    } else {
      k_list <- floor(m * c(0.05, 0.1, 0.2, 0.3, 0.5))
      p_subset <- numeric(length(k_list))
      for (i in seq_along(k_list)) {
        k <- k_list[i]
        top_k_pvals <- sort(pvals)[1:k]
        weights_k <- rep(1 / k, k)
        stat_k <- sum(weights_k * tan((0.5 - top_k_pvals) * pi))
        p_subset[i] <- 0.5 - atan(stat_k) / pi
      }
      weights2 <- rep(1 / length(p_subset), length(p_subset))
      cauchy_stat2 <- sum(weights2 * tan((0.5 - p_subset) * pi))
      p_cauchy <- 0.5 - atan(cauchy_stat2) / pi
    }
    
    print(paste0("Finished OTU: ", otu, " | p = ", signif(p_cauchy, 3)))
    
    list(pval_df = data.frame(OTU = otu, p_value = p_cauchy),
         spectrum = beta1_grid)
  }
  
  # Run without progress bar
  if (.Platform$OS.type == "unix") {
    out_list <- parallel::mclapply(otus, run_for_otu, mc.cores = parallel::detectCores())
  } else {
    out_list <- pbapply::pblapply(otus, run_for_otu)
  }
  
  # Filter valid results
  out_list <- Filter(Negate(is.null), out_list)
  results_list <- lapply(out_list, `[[`, "pval_df")
  spectra_matrix <- setNames(lapply(out_list, `[[`, "spectrum"), vapply(results_list, function(x) x$OTU, ""))
  
  result_df <- do.call(rbind, results_list)
  result_df$adj_p_value <- p.adjust(result_df$p_value, method = "bonferroni")
  
  if (!cluster) return(result_df)
  
  sig_otus <- result_df$OTU[result_df$adj_p_value < 0.05]
  if (length(sig_otus) < 2) {
    warning("Too few significant OTUs to perform clustering.")
    result_df$Cluster <- NA
    result_df$IsMedoid <- FALSE
    return(result_df)
  }
  
  spectrum_mat <- do.call(cbind, spectra_matrix[sig_otus])
  colnames(spectrum_mat) <- sig_otus
  
  dist_mat <- compute_l1_dist_matrix(spectrum_mat)
  colnames(dist_mat) <- rownames(dist_mat) <- sig_otus
  dist_obj <- as.dist(dist_mat)
  
  pamk_result <- pamk(as.matrix(dist_obj), krange = 2:min(10, length(sig_otus)), usepam = TRUE)
  k <- pamk_result$nc
  pam_result <- pam(dist_obj, k = k, diss = TRUE)
  cluster_labels <- pam_result$clustering
  medoids <- pam_result$medoids
  
  result_df$Cluster <- NA
  result_df$IsMedoid <- FALSE
  result_df$Cluster[result_df$OTU %in% names(cluster_labels)] <- cluster_labels[result_df$OTU[result_df$OTU %in% names(cluster_labels)]]
  result_df$IsMedoid[result_df$OTU %in% medoids] <- TRUE
  
  return(result_df)
}
