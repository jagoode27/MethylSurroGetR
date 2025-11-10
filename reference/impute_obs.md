# Impute Missing Observations in methyl_surro Object

This function imputes missing values in the methylation matrix of a
`methyl_surro` object using either mean or median imputation. It applies
row-wise imputation based on a specified non-missing data proportion
threshold.

## Usage

``` r
impute_obs(
  methyl_surro,
  method = c("mean", "median"),
  min_nonmiss_prop = 0,
  return_stats = FALSE,
  verbose = FALSE
)
```

## Arguments

- methyl_surro:

  An object of class `methyl_surro`, containing a methylation matrix.

- method:

  A character string indicating the imputation method. Must be either
  `"mean"` or `"median"`.

- min_nonmiss_prop:

  The minimum proportion of non-missing data required in a probe (row)
  for the imputation to proceed. Must be a numeric value between 0 and
  1.

- return_stats:

  Logical. If `TRUE`, detailed imputation statistics are added to the
  returned object. Default is `FALSE`.

- verbose:

  Logical. If `TRUE`, detailed progress messages are displayed. Default
  is `FALSE`.

## Value

A `methyl_surro` object with missing observations imputed, updating the
methylation matrix as per the given method and threshold. If
`return_stats = TRUE`, an additional `imputation_stats` component is
added containing detailed statistics about the imputation process.

## Details

The function uses vectorized operations for efficient processing of
large matrices. Probes are only imputed if they meet the minimum
non-missing data threshold and are not completely missing.

## Examples

``` r
# Load the sample data
data(beta_matrix_miss)
data(wts_df)

# Create weight vector and methyl_surro object
wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
methyl_surro_miss <- surro_set(beta_matrix_miss, wts_vec_lin, "Intercept")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.

# Apply mean imputation with a specified threshold
result <- impute_obs(methyl_surro_miss, "mean", min_nonmiss_prop = 0.5)
#> Imputed 1 values in 1 probes using mean method.
#> 1 probes were not imputed because they did not meet the threshold. Use methyl_miss() to see missing probes.

# Check the imputed result
print(result$methyl)
#>          samp1     samp2     samp3     samp4      samp5
#> cg02 0.2875775        NA        NA 0.8830174         NA
#> cg07 0.8998250 0.2460877 0.6070843 0.3279207 0.95450365
#> cg08 0.8895393 0.6928034 0.6405068 0.9942698 0.65570580
#> cg13 0.9630242 0.9022990 0.6907053 0.7954674 0.02461368
#> cg17 0.1428000 0.4145463 0.4137243 0.3688455 0.15244475
#> cg03        NA        NA        NA        NA         NA
#> cg06        NA        NA        NA        NA         NA
#> cg11        NA        NA        NA        NA         NA
#> cg15        NA        NA        NA        NA         NA
#> cg18        NA        NA        NA        NA         NA

# Get detailed imputation statistics
result_with_stats <- impute_obs(methyl_surro_miss, "mean",
                                min_nonmiss_prop = 0.5,
                                return_stats = TRUE,
                                verbose = TRUE)
#> Creating copy of methyl_surro object
#> Starting imputation analysis...
#> Found 29 missing values across 10 probes (58.0% of total values missing)
#> Using mean imputation with minimum 50.0% non-missing data threshold
#> Probe analysis:
#> - 3 probes are complete (30.0%)
#> - 5 probes are completely missing (50.0%)
#> - 2 probes have partial missing data (20.0%)
#> - 1 of 2 partial probes meet the 50.0% threshold for imputation
#> - 1 probes do not meet threshold
#> Calculating row-wise means for imputation...
#> Performing imputation...
#> Imputation completed successfully:
#> - Imputed 1 missing values across 1 probes
#> - Average imputation per probe: 1.0 values
#> - 1 probes (10.0%) were not imputed due to insufficient data
#> - Remaining missing values: 28 (56.0% of total)
#> Recommendation: 1 probes were skipped (average 40.0% completeness). Consider:
#> - Using reference_fill() to handle missing probes
#> - Lowering min_nonmiss_prop threshold to 0.3
#> - Checking data quality with methyl_miss()
#> Note: 5 completely missing probes detected. Use reference_fill() to address these before calculating surrogate values.
print(result_with_stats$imputation_stats)
#> Imputation Statistics
#> =====================
#> Method: mean
#> Threshold: 50.0% non-missing data required
#> Date: 2025-11-10 20:46:09.732577
#> 
#> Probe Summary:
#> - Total probes: 10
#> - Complete probes: 3 (30.0%)
#> - Partially missing: 2 (20.0%)
#> - Completely missing: 5 (50.0%)
#> 
#> Imputation Results:
#> - Probes imputed: 1 of 2 eligible (50.0%)
#> - Values imputed: 1 of 29 missing (3.4%)
#> - Remaining missing: 28 values
#> 
#> Skipped probes: 1 (insufficient data)
```
