# WasserSpectrum

**WasserSpectrum** is an R toolkit for distribution-aware analysis of alpha diversity. It provides:

1. **Permutation Wasserstein tests** for group comparisons  
2. **Quantile-resolved spectrum regression** to reveal where distributions differ  
3. **Post-inference Wald-type tools** (global, interval-wise, and shape contrasts)

---

## Installation

```r
# If devtools is not installed:
install.packages("devtools")

devtools::install_github("bioscinema/WasserSpectrum")
```

---

## Quick start

```r
library(WasserSpectrum)
```

### Two-group permutation Wasserstein test

```r
res2 <- wasserstein_test(
  df = df, feature_col = "Shannon", group_col = "Group",
  nperm = 2000, seed = 2025, plot = FALSE
)
res2$p_value
```

### All pairwise two-group tests across multi-level groups

```r
pair_res <- wasserstein_test_pairwise(
  df, feature_col = "Shannon", group_col = "Group",
  nperm = 2000, seed = 2025
)
subset(pair_res, p_value < 0.05)
```

### Multi-group omnibus test (Fréchet–Wasserstein ANOVA analogue)

```r
fw <- frechet_wasserstein_test(
  df, feature_col = "Shannon", group_col = "Group",
  nperm = 2000, seed = 2025, plot = FALSE
)
fw$p_value
```

### Spectrum regression (binary/numeric exposure with optional confounders)

```r
fit_bin <- wasserstein_spectrum(
  df = df, feature_col = "Shannon", outcome_col = "BMI",
  confounder_cols = c("Age", "Sex"), plot = TRUE, seed = 2025
)
fit_bin$lollipop_plot
```

### Spectrum regression (multiclass exposure with reference group)

```r
fit_mc <- wasserstein_spectrum_multiclass(
  df = df, feature_col = "Shannon",
  outcome_col = "DiseaseGroup", reference_level = "Healthy",
  confounder_cols = c("Age", "Sex"), plot = TRUE, seed = 2025
)
fit_mc$plot
```

### Post-inference: interval-wise global test over [0.2, 0.8]

```r
global <- if_manova(fit_mc, a = 0.2, b = 0.8)
global$p_value
```

### Post-inference: specific pairwise contrast over [0.2, 0.8] (Group1 vs Group2)

```r
pw <- if_manova_contrast(fit_mc, a = 0.2, b = 0.8, contrast = c(1, -1, 0))
pw$p_value
```

### Post-inference: single-curve Functional Linear Test over [0.25, 0.75]

```r
flt <- quantile_FLT(fit_bin, a = 0.25, b = 0.75)
flt$p
```

### Post-inference: shape contrast (e.g., U-shape difference between two groups)

```r
t_vec <- seq(0.1, 0.9, length.out = 5)
contrast_t <- c(1, -0.5, -1, -0.5, 1) # sums ~ 0 for pure shape
shape <- shape_contrast_test(
  spectrum_obj = fit_mc,
  t_vec = t_vec, contrast_t = contrast_t,
  contrast_group = c(1, -1, 0) # group1 vs group2
)
shape$p_value
```

---

## I. Wasserstein test for group comparison

**Goal:** Compare entire distributions (location/scale/shape) using 1D Wasserstein distance and barycenter geometry — robust to ties and skew common in microbiome data.

### `wasserstein_test`

Two-sample, permutation-based 1D Wasserstein test (Earth Mover’s Distance).

- **Inputs:** `df`, `feature_col`, `group_col` (exactly two levels), `nperm`, `seed`, `plot`  
- **Outputs:** `statistic`, `p_value`, `null_distribution`  
- **Notes:** right-tailed p-value computed as `mean(T_null >= T_obs)`.

```r
res2 <- wasserstein_test(df, "Shannon", "Group", nperm = 1000, seed = 1, plot = FALSE)
res2$statistic
res2$p_value
```

### `wasserstein_test_pairwise`

All pairwise two-sample tests across a multi-level `group_col`. Internally reuses `wasserstein_test()` on each pair.

- **Inputs:** `df`, `feature_col`, `group_col` (≥ 2 levels), `nperm`, `seed`, `pair_sep`  
- **Output:** `data.frame` with `comparison`, `statistic`, `p_value`

```r
pair_res <- wasserstein_test_pairwise(df, "Shannon", "Group", nperm = 2000, seed = 7)
pair_res[order(pair_res$p_value), ]
```

### `frechet_wasserstein_test`

Omnibus multi-group test (G > 2) using Fréchet variance in Wasserstein space (ANOVA analogue).  
Idea: compute group quantile functions on a grid, their Wasserstein barycenter, and the size-weighted sum of squared distances to that barycenter; generate a permutation null by shuffling group labels.

- **Inputs:** `df`, `feature_col`, `group_col`, `nperm`, `seed`, `plot`  
- **Outputs:** `T_obs`, `p_value`, `T_null`

```r
fw <- frechet_wasserstein_test(df, "Shannon", "Group", nperm = 2000, seed = 42, plot = FALSE)
fw$T_obs
fw$p_value
```

> **Naming note:** The omnibus function is `frechet_wasserstein_test`. Pairwise testing across groups is provided by `wasserstein_test_pairwise`.

---

## II. Spectrum analysis

**Goal:** Quantify where distributions differ along the quantile spectrum via GLS spectrum regression with B-spline basis expansion and robust (HC1) standard errors.

### `wasserstein_spectrum`

Binary or numeric exposure (optionally with confounders).

- **Inputs:**  
  `df`, `feature_col`, `outcome_col` (binary or numeric), `confounder_cols = NULL`,  
  `basis_df = 6`, `t_grid = seq(0.01, 0.99, length.out = 100)`, `alpha = 0.05`, `plot = TRUE`, `seed`  
- **Outputs:**  
  `quantiles`, `beta1`, `se1`, `lower`, `upper`, `basis`, `coef`, `vcov`, `idx_x`, `basis_df`,  
  plus `lollipop_plot` and `circular_plot` for visualization.

```r
fit_bin <- wasserstein_spectrum(
  df, feature_col = "Shannon", outcome_col = "BMI",
  confounder_cols = c("Age", "Sex"),
  basis_df = 6, alpha = 0.05, plot = TRUE, seed = 11
)
plot(fit_bin$quantiles, fit_bin$beta1, type = "l")
fit_bin$lollipop_plot
```

### `wasserstein_spectrum_multiclass`

Categorical exposure with K ≥ 2 levels; one level is set as reference. Returns group-vs-reference contrast curves βk(t) with CIs.

- **Inputs:**  
  `df`, `feature_col`, `outcome_col`, `reference_level`, `confounder_cols = NULL`,  
  `basis_df = 6`, `t_grid`, `alpha`, `plot = TRUE`, `seed`  
- **Outputs:**  
  `quantiles`, `beta` (groups × quantiles), `se`, `lower`, `upper`,  
  `comparison_levels`, `reference_level`, `basis`, `coef`, `vcov`, and `plot` (faceted)

```r
fit_mc <- wasserstein_spectrum_multiclass(
  df, "Shannon", "DiseaseGroup",
  reference_level = "Healthy",
  confounder_cols = c("Age", "Sex"),
  basis_df = 6, plot = TRUE, seed = 29
)
fit_mc$comparison_levels
fit_mc$plot
```

---

## III. Post-inference tools

Wald-type tests on integrated or shape-based features of the estimated spectrum curves.

### `if_manova`

Global interval-wise MANOVA over [a,b]: tests whether all group-vs-reference integrated effects are zero.

- **Inputs:** `spectrum_obj = wasserstein_spectrum_multiclass(...)`, `a`, `b`  
- **Outputs:** `statistic`, `df`, `p_value`, `delta`, `cov_delta`

```r
global <- if_manova(fit_mc, a = 0.2, b = 0.8)
global$statistic
global$p_value
```

### `if_manova_contrast`

User-specified contrast(s) over [a,b] (pairwise or composite).

- **Inputs:** `spectrum_obj`, `a`, `b`, `contrast` (vector or matrix; rows sum to 0 recommended)  
- **Outputs:** `statistic`, `df`, `p_value`, `contrast`, `estimate`, `cov_estimate`

```r
# Example: first two non-reference groups: 1 vs 2
pw <- if_manova_contrast(fit_mc, a = 0.2, b = 0.8, contrast = c(1, -1, 0))
pw$p_value
```

### `quantile_FLT`

Functional Linear Test for a single curve (e.g., β1(t) from `wasserstein_spectrum`) integrated over [a,b].

- **Inputs:** `spectrum_obj = wasserstein_spectrum(...)`, `a`, `b`  
- **Outputs:** `delta` (integral), `se`, `z`, `p` (two-sided)

```r
flt <- quantile_FLT(fit_bin, a = 0.25, b = 0.75)
with(flt, c(delta = delta, se = se, z = z, p = p))
```

### `shape_contrast_test`

Shape-aware testing using quantile weights (`contrast_t`) and group contrasts — ideal for U-shape vs shift comparisons.

- **Inputs:**  
  `spectrum_obj = wasserstein_spectrum_multiclass(...)`,  
  `t_vec`, `contrast_t` (should sum to ~0 for pure shape),  
  `contrast_group` (optional; rows sum to 0 recommended)  
- **Outputs:**  
  `statistic`, `df`, `p_value`, `delta_shape`, `cov_shape`,  
  `contrast_group`, `estimate`, `cov_estimate`

```r
t_vec <- seq(0.1, 0.9, length.out = 5)
contrast_t <- c(1, -0.5, -1, -0.5, 1)
shape <- shape_contrast_test(
  spectrum_obj = fit_mc,
  t_vec = t_vec, contrast_t = contrast_t,
  contrast_group = c(1, -1, 0)
)
shape$statistic
shape$p_value
```

---

