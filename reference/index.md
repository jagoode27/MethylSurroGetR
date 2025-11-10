# Package index

## Setting Up Surrogates

Initialize your surrogate object

- [`surro_set()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md)
  : Create a methyl_surro Object

## Missing Data Handling

Assess and address missing methylation data

- [`methyl_miss()`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
  : Summarize Missing Data in methyl_surro Object
- [`reference_fill()`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
  : Fill Missing Probes in methyl_surro Object Using Reference Data
- [`impute_obs()`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
  : Impute Missing Observations in methyl_surro Object

## Calculate Predictions

Compute surrogate predictions from methylation data

- [`surro_calc()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md)
  : Calculate Predicted Values from methyl_surro Object

## Data Conversion

Convert between beta values and M-values

- [`convert_m_to_beta()`](https://jagoode27.github.io/MethylSurroGetR/reference/methylConversion.md)
  [`convert_beta_to_m()`](https://jagoode27.github.io/MethylSurroGetR/reference/methylConversion.md)
  : Conversion Functions for M-Values and Beta Values

## Print Methods

S3 print methods for package objects

- [`print(`*`<imputation_stats>`*`)`](https://jagoode27.github.io/MethylSurroGetR/reference/print.imputation_stats.md)
  : Print Method for imputation_stats Objects
- [`print(`*`<methyl_miss>`*`)`](https://jagoode27.github.io/MethylSurroGetR/reference/print.methyl_miss.md)
  : Print Method for methyl_miss Objects
- [`print(`*`<reference_fill_stats>`*`)`](https://jagoode27.github.io/MethylSurroGetR/reference/print.reference_fill_stats.md)
  : Print Method for reference_fill_stats Objects

## Utility Functions

Helper functions for data generation

- [`sample_data`](https://jagoode27.github.io/MethylSurroGetR/reference/sample_data.md)
  : Sample Methylation Data for MethylSurroGetR Package

## Example Data

Sample datasets included with the package

- [`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md)
  : Complete Beta Matrix without Missing Values
- [`beta_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_miss.md)
  : Beta Matrix with Missing Values
- [`mval_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_comp.md)
  : Complete M-Values Matrix
- [`mval_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_miss.md)
  : M-Values Matrix with Missing Values
- [`wts_df`](https://jagoode27.github.io/MethylSurroGetR/reference/wts_df.md)
  : Regression Weights Data Frame
- [`ref_df`](https://jagoode27.github.io/MethylSurroGetR/reference/ref_df.md)
  : Reference Values Data Frame
