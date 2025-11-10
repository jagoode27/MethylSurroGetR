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
  calculation details

## Details

The function performs validation and processing steps:

- Validates probe-weight alignment between methylation data and weights

- Sets completely missing probes to zero

- Removes samples with any missing values

- Calculates weighted sums and applies transformations

- Provides detailed feedback about data processing

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
