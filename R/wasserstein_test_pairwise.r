#' Pairwise Permutation Wasserstein Tests (wrapper around `wasserstein_test`)
#'
#' Runs `wasserstein_test()` for all pairwise comparisons of the groups in `group_col`.
#'
#' @param df Data frame with the feature and group columns.
#' @param feature_col String; column name for the feature (e.g., diversity index).
#' @param group_col String; column name for the grouping variable (>= 2 unique values).
#' @param nperm Integer; number of permutations per pair (default 1000).
#' @param seed Integer; base seed (each pair uses `seed + i` for reproducibility).
#' @param pair_sep String; separator used in the pair label (default "-").
#'
#' @return A data.frame with columns:
#' \item{comparison}{Pair label like "group1-group2".}
#' \item{statistic}{Observed 1D Wasserstein distance for the pair.}
#' \item{p_value}{Right-tailed permutation p-value.}
#' 
#' @examples
#' \dontrun{
#' pair_res <- wasserstein_test_pairwise(df, "Shannon", "Group", nperm = 2000)
#' pair_res[order(pair_res$p_value), ]
#' }
#' @export
wasserstein_test_pairwise <- function(df, feature_col, group_col,
                                           nperm = 1000, seed = 123, pair_sep = "-") {
  # basic checks
  if (!feature_col %in% names(df)) stop("`feature_col` not found in df.")
  if (!group_col %in% names(df)) stop("`group_col` not found in df.")
  
  # drop rows with NA in required columns
  d0 <- stats::na.omit(df[, c(feature_col, group_col)])
  
  groups <- unique(d0[[group_col]])
  if (length(groups) < 2) stop("`group_col` must have at least 2 unique values.")
  
  # all pairwise combinations
  pairs <- utils::combn(groups, 2, simplify = FALSE)
  
  out <- lapply(seq_along(pairs), function(i) {
    g1 <- pairs[[i]][1]
    g2 <- pairs[[i]][2]
    
    # subset to this pair
    d_pair <- d0[d0[[group_col]] %in% c(g1, g2), , drop = FALSE]
    
    # ensure group order is consistent with label
    d_pair[[group_col]] <- droplevels(d_pair[[group_col]])
    
    # run the original test (no plotting here)
    set.seed(seed + i)
    res <- wasserstein_test(
      df = d_pair,
      feature_col = feature_col,
      group_col   = group_col,
      nperm = nperm,
      plot  = FALSE,
      seed  = seed + i
    )
    
    data.frame(
      comparison = paste0(as.character(g1), pair_sep, as.character(g2)),
      statistic  = unname(res$statistic),
      p_value    = unname(res$p_value),
      stringsAsFactors = FALSE
    )
  })
  
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
