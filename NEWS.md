# MethylSurroGetR (development version)

## New Features

* Initial package development
* Core workflow functions:
  - `surro_set()` - Create surrogate objects from methylation data and weights
  - `reference_fill()` - Fill missing probes with reference values
  - `impute_obs()` - Impute missing observations using mean or median
  - `surro_calc()` - Calculate surrogate predictions with multiple transformations
* Missing data assessment with `methyl_miss()`
* Beta value and M-value conversion functions:
  - `convert_beta_to_m()` - Convert beta values to M-values
  - `convert_m_to_beta()` - Convert M-values to beta values
* Comprehensive input validation and error messages
* Detailed diagnostic reporting options
* Support for linear, probability (logistic), and count (Poisson) transformations
* Example datasets for testing and learning

## Documentation

* Added "Getting Started with MethylSurroGetR" vignette
* Complete function documentation with examples
* Package website at <https://jagoode27.github.io/MethylSurroGetR/>