# Calculate Predicted Values from methyl_surro Object

This function calculates predicted values based on the provided
methyl_surro object. If any probes (rows) are completely missing, their
values are set to zero. Samples with any missing values are removed from
the calculation. The function aligns and multiplies the weights with the
methylation matrix, adds the intercept if specified, and applies the
specified transformation to the results.

## Usage

``` r
surro_calc(
  methyl_surro,
  transform = c("linear", "count", "probability"),
  verbose = FALSE,
  return_diagnostics = FALSE
)
```

## Arguments

- methyl_surro:

  An object of class `methyl_surro`.

- transform:

  A character string specifying the transformation to apply. Can be
  `"linear"`, `"count"`, or `"probability"`.

- verbose:

  Logical. If `TRUE`, detailed progress and diagnostic messages are
  displayed. Default is `FALSE`.

- return_diagnostics:

  Logical. If `TRUE`, returns additional diagnostic information. Default
  is `FALSE`.

## Value

If `return_diagnostics = FALSE`: A named vector of predicted values with
names corresponding to the sample names. If `return_diagnostics = TRUE`:
A list containing:

- predictions:

  The named vector of predicted values

- diagnostics:

  A list of diagnostic information including missing data handling and
  calculation details (see Diagnostics section)

## Details

The function performs validation and processing steps:

- Validates probe-weight alignment between methylation data and weights

- Sets completely missing probes to zero

- Removes samples with any missing values

- Calculates weighted sums and applies transformations

- Provides detailed feedback about data processing

## Transformation Equations

The surrogate biomarker value for each sample is calculated using the
following equations, where \\X\\ is the methylation matrix, \\w\\ is the
weight vector, and \\b\\ is the intercept:

**Linear transformation (`transform = "linear"`):** \$\$y = X^T w +
b\$\$ This is the raw weighted sum of methylation values, suitable for
continuous outcomes such as estimated cell counts or continuous
biomarker values.

**Count transformation (`transform = "count"`):** \$\$y = \exp(X^T w +
b)\$\$ This applies an exponential transformation to the linear
predictor, appropriate for count-based outcomes following a log-linear
model (e.g., Poisson regression results). Values are constrained to be
positive.

**Probability transformation (`transform = "probability"`):** \$\$y =
\frac{1}{1 + \exp(-(X^T w + b))}\$\$ This applies the inverse logit
(logistic) transformation, converting the linear predictor to a
probability scale (range: 0 to 1). Appropriate for binary outcomes or
risk scores derived from logistic regression models.

## Diagnostics

When `return_diagnostics = TRUE`, the function returns a list where the
`diagnostics` component contains detailed information about the
calculation:

- transform:

  Character string indicating which transformation was applied

- n_samples_original:

  Number of samples in the input methylation matrix

- n_probes_original:

  Number of probes (rows) in the input methylation matrix

- n_weights:

  Number of weights provided

- intercept_used:

  Logical indicating whether an intercept was included

- n_common_probes:

  Number of probes present in both weights and methylation data

- n_missing_in_weights:

  Number of probes in methylation data without corresponding weights

- n_missing_in_methyl:

  Number of weight probes not found in methylation data

- n_completely_missing_probes:

  Number of probes with all NA values (set to zero)

- n_samples_with_missing:

  Number of samples containing any NA values (removed)

- probes_set_to_zero:

  Same as n_completely_missing_probes

- samples_removed:

  Same as n_samples_with_missing

- n_final_samples:

  Number of samples in the final prediction vector

- final_result_range:

  Named numeric vector with min and max values of predictions

## Examples

``` r
# Load sample data
data(beta_matrix_comp, package = "MethylSurroGetR")
data(wts_df, package = "MethylSurroGetR")

# Create weight vector and methyl_surro object
wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
methyl_surro_comp <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.

# Basic calculation
pred_lin <- surro_calc(methyl_surro_comp, transform = "linear")
#> 5 probes set to zero. Use reference_fill().
print(pred_lin)
#>    samp1    samp2    samp3    samp4    samp5 
#> 1.196622 1.199377 1.205871 1.198700 1.197060 

# With diagnostic information
result_detailed <- surro_calc(methyl_surro_comp,
                              transform = "linear",
                              return_diagnostics = TRUE,
                              verbose = TRUE)
#> Starting surrogate calculation with linear transformation...
#> Input: 10 probes x 5 samples
#> Weights: 10 probe weights available
#> Intercept: 1.211000
#> 5 probes set to zero. Use reference_fill().
#> Calculating weighted sum...
#> Added intercept: 1.211000
#> Applying linear transformation...
#> Calculation completed successfully:
#> - Final sample count: 5
#> - Result range: 1.1966 to 1.2059
print(result_detailed$diagnostics)
#> $transform
#> [1] "linear"
#> 
#> $n_samples_original
#> [1] 5
#> 
#> $n_probes_original
#> [1] 10
#> 
#> $n_weights
#> [1] 10
#> 
#> $intercept_used
#> [1] TRUE
#> 
#> $n_common_probes
#> [1] 10
#> 
#> $n_missing_in_weights
#> [1] 0
#> 
#> $n_missing_in_methyl
#> [1] 0
#> 
#> $n_completely_missing_probes
#> [1] 5
#> 
#> $n_samples_with_missing
#> [1] 0
#> 
#> $probes_set_to_zero
#> [1] 5
#> 
#> $samples_removed
#> [1] 0
#> 
#> $n_final_samples
#> [1] 5
#> 
#> $final_result_range
#>      min      max 
#> 1.196622 1.205871 
#> 
```
