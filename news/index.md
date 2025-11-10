# Changelog

## MethylSurroGetR 0.0.0.9000

### New Features

- Initial package development
- Core workflow functions:
  - [`surro_set()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md) -
    Create surrogate objects from methylation data and weights
  - [`reference_fill()`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md) -
    Fill missing probes with reference values
  - [`impute_obs()`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md) -
    Impute missing observations using mean or median
  - [`surro_calc()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md) -
    Calculate surrogate predictions with multiple transformations
- Missing data assessment with
  [`methyl_miss()`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
- Beta value and M-value conversion functions:
  - [`convert_beta_to_m()`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_beta_to_m.md) -
    Convert beta values to M-values
  - [`convert_m_to_beta()`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_m_to_beta.md) -
    Convert M-values to beta values
- Comprehensive input validation and error messages
- Detailed diagnostic reporting options
- Support for linear, probability (logistic), and count (Poisson)
  transformations
- Example datasets for testing and learning

### Documentation

- Added “Getting Started with MethylSurroGetR” vignette
- Complete function documentation with examples
- Package website at <https://jagoode27.github.io/MethylSurroGetR/>
