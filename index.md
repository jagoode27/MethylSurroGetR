# MethylSurroGetR

### 🚧 Under Development 🚧

**This package is still in development and not yet ready for general
use. Proceed with caution!**

The goal of MethylSurroGetR is to offer user-friendly tools and
functions for generating predicted values from existing DNA methylation
surrogates. It aims to enhance the accessibility and reproducibility of
data analysis with methylation surrogates.

## Installation

You can install the development version of MethylSurroGetR from GitHub
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
wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))

# Create surrogate object
my_surro <- surro_set(
  methyl = beta_matrix_comp,
  weights = wts_vec,
  intercept = "Intercept"
)

# Fill missing probes with reference values
ref_vec <- setNames(ref_df$mean, rownames(ref_df))
my_surro <- reference_fill(my_surro, reference = ref_vec)

# Calculate predictions
predictions <- surro_calc(my_surro, transform = "linear")
```

## Getting Help

See the [Getting Started
vignette](https://jagoode27.github.io/MethylSurroGetR/articles/MethylSurroGetR.html)
for detailed examples and workflows.

For bug reports and feature requests, please [open an
issue](https://github.com/jagoode27/MethylSurroGetR/issues).
