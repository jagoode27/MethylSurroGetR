# MethylSurroGetR

MethylSurroGetR is an R package designed to facilitate the construction
of DNA methylation surrogate biomarker values from existing studies. It
provides comprehensive tools for handling missing methylation data
through various imputation strategies, including mean/median imputation
for missing observations and reference-based filling for missing probes.
The package supports flexible surrogate calculation with “linear”,
“probability” (logistic), and “count” (Poisson) transformations,
comprehensive input validation, and detailed diagnostic reporting.
Additional utilities include conversion between beta values and
M-values, as well as missing data analysis.

## Installation

You can install the most current version of MethylSurroGetR from GitHub
with:

``` r

# install.packages("remotes")
remotes::install_github("jagoode27/MethylSurroGetR")
```

## Example

This is a basic example workflow:

``` r

library(MethylSurroGetR)

# Load example data
data(beta_matrix_comp)
data(wts_df)
data(ref_df)

# Create weights vector
wts_vec_lin <- with(wts_df, setNames(wt_lin, rownames(wts_df)))[!is.na(wts_df$wt_lin)]

# Create surrogate object
my_surro <- surro_set(
  methyl = beta_matrix_comp,
  weights = wts_vec_lin,
  intercept = "Intercept"
)

# Fill missing probes with reference values
ref_vec_mean <- with(ref_df, setNames(mean, rownames(ref_df)))[!is.na(ref_df$mean)]
my_surro <- reference_fill(my_surro, reference = ref_vec_mean)

# Calculate predictions
predictions <- surro_calc(my_surro, transform = "linear")
```

## Getting Help

See the [Getting Started
vignette](https://jagoode27.github.io/MethylSurroGetR/articles/MethylSurroGetR.html)
for detailed examples and workflows.

For bug reports and feature requests, please [open an
issue](https://github.com/jagoode27/MethylSurroGetR/issues).
