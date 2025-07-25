# WasserSpectrum 

`WasserSpectrum` provides a robust framework for quantile-wise inference of group differences in univariate outcomes, leveraging the Wasserstein distance and Generalized Least Squares (GLS). It is particularly useful for analyzing skewed, zero-inflated, or heavy-tailed data commonly encountered in microbial diversity, ecological studies, and biomedical data.

Key features:

* Permutation-based global tests using Wasserstein distances
* Spectrum regression models estimating functional effects $\beta(t)$
* Quantile-wise association testing using Cauchy combination methods (single or multiple features)
* Interval-wise Wald tests for effect integration
* Shape contrast tests for pattern comparison
* Quantile-wise multivariate hypothesis testing (MANOVA)
* Flexible plotting and visualization tools

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("bioscinema/WasserSpectrum")
library(WasserSpectrum)
```

---

## Core Functions and Usage

### 1. Wasserstein-based Permutation Tests

* `wasserstein_test()`
  Performs a **two-group** permutation test based on 1D Wasserstein distance. Robust to outliers and skewed distributions.

* `frechet_wasserstein_test()`
  Extends the above to **multi-group** settings ($G \geq 2$) using Fréchet mean and variance.

```r
wasserstein_test(df, feature_col = "Shannon", group_col = "Group")
frechet_wasserstein_test(df, feature_col = "Shannon", group_col = "Group")
```

---

### 2. Functional Spectrum Estimation (GLS)

* `wasserstein_spectrum()`
  Estimates a smooth effect curve $\beta(t)$ for **binary or numeric** exposures using logistic regression at multiple quantiles and B-spline basis.

* `wasserstein_spectrum_multiclass()`
  Extends spectrum estimation to **multi-class categorical** exposures. Supports covariate adjustment.

```r
fit <- wasserstein_spectrum(df, feature_col = "Shannon", outcome_col = "Group")
plot(fit$plot)  # Returns a ggplot object
```

---

### 2b. Quantile-Based Association Testing Using the Wasserstein Spectrum

* `wasserstein_quantile_single()`
  Performs quantile-wise inference on a single feature (e.g., OTU) using a B-spline projection of the outcome effect curve $\beta(t)$.
Aggregates quantile-wise p-values using Cauchy combination (with optional adaptive version).
Returns the global p-value, area under the effect curve (AUC), and absolute AUC.

```r
wasserstein_quantile_single(
  df = t(count_data),         # features in columns, samples in rows
  feature_col = "G000830295", # the feature to test
  metadata = metadata,
  outcome_col = "BMI",
  confounder_cols = c("Age", "Gender"),
  adaptive = FALSE
)
```

* `wasserstein_quantile_test()`
  Applies *wasserstein_quantile_single()* to each *feature (column in df)* and performs parallelized quantile-wise inference across all features.
Returns a data frame with global p-values, AUCs, and Bonferroni-adjusted p-values for multiple testing correction.

```r
wasserstein_quantile_test(
  df = t(count_data),  # function expects features in columns, samples in rows
  metadata = metadata,
  outcome_col = "BMI",
  confounder_cols = c("Age", "Gender"),
  adaptive = FALSE
)
```

---

### 3. Quantile-Level MANOVA & Shape Tests

* `spectrum_manova()`
  Computes F-statistic and p-value across quantiles for multivariate Wald tests on all $\beta_k(t)$.

* `plot_manova()`
  Visualizes the test result from `spectrum_manova()` as F-values and $-\log_{10}(p)$ across quantiles.

* `shape_contrast_test()`
  Tests **shape differences** in functional effect curves between groups. Can compare U-shaped vs. flat profiles.

```r
manova_out <- spectrum_manova(fit)
plot(manova_out)  # Plot F and p curves
```

---

### 4. Interval-Based Integration Tests

* `quantile_FLT()`
  For binary exposures, performs a **Functional Linear Test (FLT)** by integrating $\beta(t)$ over a specified interval $[a,b]$.

* `if_manova()`
  Generalizes FLT for **multi-class** outcomes. Computes Wald-type statistic for $\int_a^b \beta_k(t) dt$.

* `if_manova_contrast()`
  Extends `if_manova()` to handle **custom contrasts** among group-level integrals.

```r
if_manova(fit, a = 0.1, b = 0.9)  # Global interval test
if_manova_contrast(fit, a = 0.2, b = 0.6, contrast = c(1, -1, 0))  # Group 1 vs Group 2
```

---

## Full Workflow

1. **Run global Wasserstein test** to assess if overall group distributions differ.
2. **Fit spectrum model** with `wasserstein_spectrum[_multiclass]()`.
3. **Visualize** estimated effect curves $\beta_k(t)$ and confidence bands.
4. **Test hypotheses** on:

   * Quantile range (via `quantile_FLT`, `if_manova`, `if_manova_contrast`)
   * Shape contrasts (`shape_contrast_test`)
   * Full-grid profile (`spectrum_manova` + `plot.manova`)
   * Single/multiple feature association (`wasserstein_quantile_single`, `wasserstein_quantile_test`)
---

## License

This package is released under **GPL (>= 2)**. See LICENSE file for details.

---

## Project URL

[https://github.com/bioscinema/WasserSpectrum](https://github.com/bioscinema/WasserSpectrum)

Bug reports and contributions are welcome via GitHub Issues or Pull Requests.
