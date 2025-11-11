# Fill Missing Probes in methyl_surro Object Using Reference Data

This function fills in missing probes in the methylation matrix of a
`methyl_surro` object based on external reference data. Users can
specify if imputation should consider only fully missing rows, partially
missing rows, or the entire dataset.

## Usage

``` r
reference_fill(
  methyl_surro,
  reference,
  type = c("probes", "obs", "all"),
  return_stats = FALSE,
  verbose = FALSE
)
```

## Arguments

- methyl_surro:

  An object of class `methyl_surro`, containing a methylation matrix.

- reference:

  A named numeric vector containing reference values for probes to
  impute missing data.

- type:

  A character string specifying the scope for filling missing data:

  - `"probes"`: Fill in only missing probes (completely missing rows).

  - `"obs"`: Fill in only missing observations (partially missing
    values).

  - `"all"`: Fill in missing probes and observations.

- return_stats:

  Logical. If `TRUE`, detailed filling statistics are added to the
  returned object. Default is `FALSE`.

- verbose:

  Logical. If `TRUE`, detailed progress messages are displayed. Default
  is `FALSE`.

## Value

A `methyl_surro` object with its methylation matrix updated by filling
in the specified missing probes based on the reference data. If
`return_stats = TRUE`, an additional `reference_fill_stats` component is
added containing detailed statistics about the filling process.

## Details

The function uses vectorized operations for efficient processing of
large matrices. Reference values are validated to ensure they fall
within reasonable ranges for methylation data (0-1 for beta values,
typically -10 to 10 for M-values).

**Filling Strategy:** The function applies reference values consistently
across all samples for each probe, which is appropriate for
population-level reference data but may not capture individual-level
variation.

## Examples

``` r
# Load the sample data
data(beta_matrix_miss)
data(wts_df)
data(ref_df)

# Create weight vector and methyl_surro object
wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
methyl_surro_miss <- surro_set(beta_matrix_miss, wts_vec_lin, "Intercept")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.

# Extract specific reference columns as named vectors
ref_mean <- setNames(ref_df$mean, rownames(ref_df))
ref_median <- setNames(ref_df$median, rownames(ref_df))

# Apply reference filling using reference vector
result <- reference_fill(methyl_surro = methyl_surro_miss,
                         reference = ref_mean,
                         type = "probes",
                         verbose = TRUE)
#> Starting reference filling analysis...
#> Matrix dimensions: 10 probes x 5 samples
#> Complete probes: 3 (30.0%)
#> Partially missing probes: 2 (20.0%)
#> Completely missing probes: 5 (50.0%)
#> Reference data available for 20 probes
#> Filling 5 completely missing probes with reference values...
#> Filled 25 values using reference data (probes strategy).
#> Reference filling completed successfully:
#> - Filled 5 completely missing probes (25 values)
#> - Total values filled: 25 (86.2% of originally missing)
#> - Remaining missing values: 4 (8.0% of total)

# Check the result after filling
print(result$methyl)
#>          samp1     samp2     samp3     samp4      samp5
#> cg02 0.2875775        NA        NA 0.8830174         NA
#> cg07 0.8998250 0.2460877        NA 0.3279207 0.95450365
#> cg08 0.8895393 0.6928034 0.6405068 0.9942698 0.65570580
#> cg13 0.9630242 0.9022990 0.6907053 0.7954674 0.02461368
#> cg17 0.1428000 0.4145463 0.4137243 0.3688455 0.15244475
#> cg03 0.4948262 0.4948262 0.4948262 0.4948262 0.49482616
#> cg06 0.5526592 0.5526592 0.5526592 0.5526592 0.55265924
#> cg11 0.4566024 0.4566024 0.4566024 0.4566024 0.45660238
#> cg15 0.4004940 0.4004940 0.4004940 0.4004940 0.40049405
#> cg18 0.3923206 0.3923206 0.3923206 0.3923206 0.39232059

# Get detailed filling statistics
result_with_stats <- reference_fill(methyl_surro_miss, ref_median,
                                     type = "all",
                                     return_stats = TRUE, verbose = TRUE)
#> Starting reference filling analysis...
#> Matrix dimensions: 10 probes x 5 samples
#> Complete probes: 3 (30.0%)
#> Partially missing probes: 2 (20.0%)
#> Completely missing probes: 5 (50.0%)
#> Reference data available for 20 probes
#> Filling 5 completely missing probes with reference values...
#> Filling missing observations in 2 partially missing probes...
#> Filled 29 values using reference data (all strategy).
#> Reference filling completed successfully:
#> - Filled 5 completely missing probes (25 values)
#> - Filled missing observations in 2 probes (4 values)
#> - Total values filled: 29 (100.0% of originally missing)
#> - No missing values remaining
print(result_with_stats$reference_fill_stats)
#> Reference Filling Statistics
#> ============================
#> Type: all
#> Date: 2025-11-11 01:13:34.127848
#> 
#> Matrix Summary:
#> - Total probes: 10
#> - Total samples: 5
#> - Complete probes: 10 (100.0%)
#> - Partially missing: 0 (0.0%)
#> - Completely missing: 0 (0.0%)
#> 
#> Reference Data:
#> - Available reference probes: 20
#> - Reference value range: 0.266 to 0.795
#> 
#> Filling Results:
#> - Complete probes filled: 5 (25 values)
#> - Partial probes filled: 2 (4 values)
#> - Total values filled: 29 of 29 missing (100.0%)
#> - No missing values remaining
```
