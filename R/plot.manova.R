#' Plot Quantile-wise MANOVA Results
#'
#' Visualizes the output of \code{\link{spectrum_manova}} as a function of quantile level,
#' showing either the F-statistics, \(-\log_{10}(p)\)-values, or both stacked.
#'
#' @param manova_df A data frame returned by \code{\link{spectrum_manova}}, containing columns
#' \code{"quantile"}, \code{"Fval"}, and \code{"pval"}.
#' @param type Type of plot to return: \code{"F"} for F-statistics only,
#' \code{"p"} for \(-\log_{10}(p)\)-values only, or \code{"both"} (default) for stacked plots.
#' @param p_cutoff Significance threshold to highlight in the \(-\log_{10}(p)\) plot (default is 0.05).
#' @param line_color Color for the line plots (default is \code{"darkblue"}).
#' @param ribbon_color Ignored in current implementation; reserved for future use.
#' @param title Title for the plot (default is \code{"Quantile-wise MANOVA"}).
#' @param base_size Base font size for ggplot theme (default is 14).
#'
#' @return A \code{ggplot} object. If \code{type = "both"}, returns a vertically stacked plot using \code{patchwork}.
#'
#' @details
#' This function visualizes the test results from \code{\link{spectrum_manova}} to help identify
#' which quantile ranges show significant multigroup separation. The F-statistics reflect the
#' Hotelling-type test magnitude, while \(-\log_{10}(p)\)-values help assess significance relative
#' to a chosen cutoff (default 0.05).
#'
#' If \code{type = "both"}, the function returns a two-panel plot with F-statistics on top and
#' \(-\log_{10}(p)\)-values below.
#'
#' @examples
#' \dontrun{
#' manova_result <- spectrum_manova(fit)
#' plot.manova(manova_result, type = "both", p_cutoff = 0.01)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal geom_hline
#' @importFrom patchwork plot_layout
#' @export
plot.manova <- function(manova_df,
                        type = c("both", "F", "p"),
                        p_cutoff = 0.05,
                        line_color = "darkblue",
                        ribbon_color = "lightblue",
                        title = "Quantile-wise MANOVA",
                        base_size = 14) {
  library(ggplot2)
  type <- match.arg(type)
  
  if (!all(c("quantile", "Fval", "pval") %in% colnames(manova_df))) {
    stop("Input must contain 'quantile', 'Fval', and 'pval' columns.")
  }
  
  gF <- ggplot(manova_df, aes(x = quantile, y = Fval)) +
    geom_line(color = line_color, linewidth = 1.2) +
    labs(y = "F-statistic", x = "Quantile level t",
         title = paste(title, "(F values)")) +
    theme_minimal(base_size = base_size)
  
  gp <- ggplot(manova_df, aes(x = quantile, y = -log10(pval))) +
    geom_line(color = line_color, linewidth = 1.2) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "gray50") +
    labs(y = expression(-log[10](p(t))), x = "Quantile level t",
         title = paste(title, "(-log10 p-values)")) +
    theme_minimal(base_size = base_size)
  
  if (type == "F") return(gF)
  if (type == "p") return(gp)
  
  # type == "both"
  suppressPackageStartupMessages(library(patchwork))
  return(gF / gp + plot_layout(heights = c(1, 1)))
}