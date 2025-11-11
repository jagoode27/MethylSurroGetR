# MethylSurroGetR: DNA Methylation Surrogate Calculations

MethylSurroGetR provides a comprehensive toolkit for calculating DNA
methylation surrogate biomarkers from existing studies. The package
handles the complete workflow from data preparation through prediction
calculation, with robust handling of missing data and flexible
transformation options.

## Main Workflow Functions

The typical workflow involves four main steps:

- [`surro_set`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md):

  Create a surrogate object by combining methylation data with surrogate
  weights

- [`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md):

  Fill missing CpG probes using reference values (e.g., from population
  means)

- [`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md):

  Impute missing observations within samples using mean or median
  imputation

- [`surro_calc`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md):

  Calculate surrogate predictions with linear, logistic (probability),
  or Poisson (count) transformations

## Data Assessment Functions

- [`methyl_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md):

  Comprehensive missing data assessment for methylation matrices

## Data Conversion Functions

- [`convert_beta_to_m`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_beta_to_m.md):

  Convert beta values (0-1 scale) to M-values (log-ratio scale)

- [`convert_m_to_beta`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_m_to_beta.md):

  Convert M-values back to beta values

## Sample Datasets

The package includes example datasets for learning and testing:

- [`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md):

  Complete beta value matrix (15 probes × 5 samples)

- [`beta_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_miss.md):

  Beta matrix with missing values

- [`mval_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_comp.md):

  Complete M-value matrix

- [`mval_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_miss.md):

  M-value matrix with missing values

- [`wts_df`](https://jagoode27.github.io/MethylSurroGetR/reference/wts_df.md):

  Example surrogate weights with three transformation types

- [`ref_df`](https://jagoode27.github.io/MethylSurroGetR/reference/ref_df.md):

  Reference values for missing probe imputation

## Getting Started

To get started with MethylSurroGetR:

1.  Load your methylation data (beta or M-values) as a matrix with CpG
    probes as rows and samples as columns

2.  Obtain surrogate weights from a published study or your own model

3.  Follow the basic workflow:

          # Create surrogate object
          my_surro <- surro_set(methyl_matrix, weights_vector, intercept = "Intercept")

          # Handle missing probes
          my_surro <- reference_fill(my_surro, reference_values)

          # Impute missing observations (if needed)
          my_surro <- impute_obs(my_surro, method = "mean")

          # Calculate predictions
          predictions <- surro_calc(my_surro, transform = "linear")
          

4.  See
    [`vignette("MethylSurroGetR")`](https://jagoode27.github.io/MethylSurroGetR/articles/MethylSurroGetR.md)
    for detailed examples

## Key Features

- Flexible missing data handling with multiple strategies

- Support for linear, logistic (probability), and Poisson (count)
  transformations

- Comprehensive input validation and informative error messages

- Detailed diagnostic reporting for transparency

- Memory-efficient operations for large-scale datasets

- Extensive test coverage ensuring reliability

## Authors

Joshua A. Goode <jagoode@umich.edu> (ORCID: 0000-0003-3290-0284)

## See Also

Useful links:

- Package website: <https://jagoode27.github.io/MethylSurroGetR/>

- GitHub repository: <https://github.com/jagoode27/MethylSurroGetR>

- Report bugs: <https://github.com/jagoode27/MethylSurroGetR/issues>

## See also

Useful links:

- <https://github.com/jagoode27/MethylSurroGetR>

- <https://jagoode27.github.io/MethylSurroGetR/>

- Report bugs at <https://github.com/jagoode27/MethylSurroGetR/issues>

## Author

**Maintainer**: Joshua A. Goode <jagoode@umich.edu>
([ORCID](https://orcid.org/0000-0003-3290-0284))
