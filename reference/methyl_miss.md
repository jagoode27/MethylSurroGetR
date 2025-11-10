# Summarize Missing Data in methyl_surro Object

This function summarizes the missing data in a methyl_surro object. It
identifies probes with partial missing data (missing observations) and
probes missing from all samples (missing probes), along with
comprehensive summary statistics.

## Usage

``` r
methyl_miss(methyl_surro)
```

## Arguments

- methyl_surro:

  An object of class `methyl_surro`.

## Value

A list of class `methyl_miss` with the following elements:

- missing_obs:

  A named numeric vector with the proportion of missing observations for
  each probe where this value is greater than 0 but less than 1.

- missing_probes:

  A character vector of probes that are missing in all samples.

- summary:

  A list containing summary statistics about the missing data patterns.

The `summary` component contains:

- total_probes:

  Total number of probes in the methylation matrix.

- total_samples:

  Total number of samples in the methylation matrix.

- n_complete_probes:

  Number of probes with no missing values.

- n_missing_obs:

  Number of probes with partial missing data.

- n_missing_probes:

  Number of completely missing probes.

- overall_missing_rate:

  Proportion of all matrix values that are missing.

- missing_obs_rate:

  Proportion of probes that have partial missing data.

- missing_probes_rate:

  Proportion of probes that are completely missing.

- complete_probes_rate:

  Proportion of probes that are complete.

## Examples

``` r
# Load Methylation Beta Matrix
data(beta_matrix_miss, package = "MethylSurroGetR")

# Load Weights Data Frame and create vector
data(wts_df, package = "MethylSurroGetR")
wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))

# Build the methyl_surro Object
surrogate <- surro_set(methyl = beta_matrix_miss,
                       weights = wts_vec_lin,
                       intercept = "Intercept")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.

# Summarizing Missing Values
missing_summary <- methyl_miss(methyl_surro = surrogate)
print(missing_summary)
#> Missing Data Summary for methyl_surro Object
#> ============================================
#> 
#> Total probes: 10
#> Total samples: 5
#> Complete probes: 3 (30.0%)
#> Probes with missing observations: 2 (20.0%)
#> Completely missing probes: 5 (50.0%)
#> Overall missing rate: 58.0%
#> 
#> Probes with partial missing data:
#> cg02 cg07 
#>  0.6  0.2 
#> 
#> Completely missing probes:
#> cg03, cg06, cg11, cg15, cg18

# Access specific components
missing_summary$missing_obs
#> cg02 cg07 
#>  0.6  0.2 
missing_summary$summary$overall_missing_rate
#> [1] 0.58
```
