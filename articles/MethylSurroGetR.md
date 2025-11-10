# Getting Started with MethylSurroGetR

## Introduction

DNA methylation surrogates are statistical models that predict
biological or clinical outcomes using methylation values from specific
CpG sites. **MethylSurroGetR** provides a comprehensive toolkit for
working with these surrogates, handling common challenges like missing
data and providing flexible calculation methods.

This vignette demonstrates the main features of the package using
example data.

### Installation

``` r
library(MethylSurroGetR)
```

## Basic Workflow

The typical workflow involves four main steps:

1.  **Create a surrogate object** with
    [`surro_set()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md)
2.  **Handle missing probes** with
    [`reference_fill()`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
3.  **Impute missing observations** with
    [`impute_obs()`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
    (if needed)
4.  **Calculate predictions** with
    [`surro_calc()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md)

Let’s walk through each step with examples.

## Working with Complete Data

### Step 1: Load Your Data

The package includes sample datasets for demonstration:

``` r
# Load example methylation data (beta values)
data(beta_matrix_comp)

# Load surrogate weights
data(wts_df)

# Load reference values for missing probes
data(ref_df)

# Preview the data
head(beta_matrix_comp)
#>          samp1     samp2      samp3     samp4     samp5
#> cg01 0.1028646 0.3203732 0.48290240 0.7205963 0.3694889
#> cg02 0.2875775 0.7883051 0.40897692 0.8830174 0.9404673
#> cg04 0.4348927 0.1876911 0.89035022 0.1422943 0.9842192
#> cg05 0.9849570 0.7822943 0.91443819 0.5492847 0.1542023
#> cg07 0.8998250 0.2460877 0.04205953 0.3279207 0.9545036
#> cg08 0.8895393 0.6928034 0.64050681 0.9942698 0.6557058
```

The methylation matrix should have: - CpG sites as **rows** (with CpG
names as row names) - Samples as **columns** - Beta values (0-1 scale)
as entries

### Step 2: Create a Surrogate Object

Use
[`surro_set()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md)
to combine your methylation data with surrogate weights:

``` r
# Extract weights as a named vector
wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))

# Create surrogate object
my_surro <- surro_set(
  methyl = beta_matrix_comp,
  weights = wts_vec,
  intercept = "Intercept"
)
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.

# View the structure
str(my_surro)
#> List of 3
#>  $ methyl   : num [1:10, 1:5] 0.288 0.9 0.89 0.963 0.143 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:10] "cg02" "cg07" "cg08" "cg13" ...
#>   .. ..$ : chr [1:5] "samp1" "samp2" "samp3" "samp4" ...
#>  $ weights  : Named num [1:10] -0.00908 -0.00116 0.00598 -0.00756 0.00122 ...
#>   ..- attr(*, "names")= chr [1:10] "cg02" "cg03" "cg06" "cg07" ...
#>  $ intercept: Named num 1.21
#>   ..- attr(*, "names")= chr "Intercept"
#>  - attr(*, "class")= chr "methyl_surro"
```

**Key points:** - `weights` should be a named numeric vector with CpG
names - `intercept` specifies which weight element is the intercept
term - The function automatically aligns probes and adds missing ones as
NA

### Step 3: Fill Missing Probes

If your data is missing probes needed by the surrogate, use reference
values:

``` r
# Extract reference values
ref_vec <- setNames(ref_df$mean, rownames(ref_df))

# Fill missing probes with reference values
my_surro <- reference_fill(
  methyl_surro = my_surro,
  reference = ref_vec,
  type = "probes"  # Only fill completely missing probes
)
#> Filled 25 values using reference data (probes strategy).
```

**Fill types:**

- `"probes"` - Only fill completely missing probes (recommended)
- `"obs"` - Only fill individual missing observations
- `"all"` - Fill both missing probes and observations

### Step 4: Calculate Predictions

Now calculate the surrogate predictions:

``` r
# Calculate linear predictions
predictions <- surro_calc(
  methyl_surro = my_surro,
  transform = "linear"
)

print(predictions)
#>    samp1    samp2    samp3    samp4    samp5 
#> 1.197718 1.200473 1.206967 1.199796 1.198156
```

**Available transformations:** - `"linear"` - No transformation
(default) - `"probability"` - Logistic transformation (1 / (1 +
exp(-x))) - `"count"` - Exponential transformation (exp(x))

## Handling Missing Data

### Checking for Missing Data

Use
[`methyl_miss()`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
to assess missing data patterns:

``` r
# Load data with missing values
data(beta_matrix_miss)

# Create surrogate object
surro_miss <- surro_set(
  methyl = beta_matrix_miss,
  weights = wts_vec,
  intercept = "Intercept"
)
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.

# Analyze missing data
missing_info <- methyl_miss(surro_miss)
print(missing_info)
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
```

This shows:

- Number of complete probes
- Probes with partial missingness
- Completely missing probes
- Overall missing data rate

### Imputing Missing Observations

For probes with partial missingness, use
[`impute_obs()`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md):

``` r
# Impute using row means
surro_imputed <- impute_obs(
  methyl_surro = surro_miss,
  method = "mean",
  min_nonmiss_prop = 0.5,  # Only impute if ≥50% observed
  verbose = TRUE
)
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
```

**Parameters:**

- `method` - “mean” or “median” imputation
- `min_nonmiss_prop` - Minimum proportion of non-missing values required
  (0-1)
- `return_stats` - Get detailed imputation statistics

### Complete Workflow with Missing Data

``` r
# Start with missing data
data(beta_matrix_miss)
data(wts_df)
data(ref_df)

# Create vectors
wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
ref_vec <- setNames(ref_df$mean, rownames(ref_df))

# Full pipeline
predictions_miss <- beta_matrix_miss |>
  surro_set(weights = wts_vec, intercept = "Intercept") |>
  reference_fill(reference = ref_vec, type = "probes") |>
  impute_obs(method = "mean", min_nonmiss_prop = 0.5) |>
  surro_calc(transform = "linear")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.
#> Filled 25 values using reference data (probes strategy).
#> Imputed 1 values in 1 probes using mean method.
#> 1 probes were not imputed because they did not meet the threshold. Use methyl_miss() to see missing probes.
#> 3 samples omitted. Use impute_obs().

print(predictions_miss)
#>    samp1    samp4 
#> 1.197718 1.199796
```

## Different Surrogate Types

### Linear Surrogates

For continuous outcomes (e.g., estimated white blood cell counts):

``` r
wts_linear <- setNames(wts_df$wt_lin, rownames(wts_df))

predictions_linear <- beta_matrix_comp |>
  surro_set(weights = wts_linear, intercept = "Intercept") |>
  reference_fill(reference = ref_vec, type = "probes") |>
  surro_calc(transform = "linear")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.
#> Filled 25 values using reference data (probes strategy).

print(predictions_linear)
#>    samp1    samp2    samp3    samp4    samp5 
#> 1.197718 1.200473 1.206967 1.199796 1.198156
```

### Probability Surrogates

For binary outcomes (e.g., disease risk):

``` r
wts_prob <- setNames(wts_df$wt_prb, rownames(wts_df))

predictions_prob <- beta_matrix_comp |>
  surro_set(weights = wts_prob, intercept = "Intercept") |>
  reference_fill(reference = ref_vec, type = "probes") |>
  surro_calc(transform = "probability")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.
#> Filled 25 values using reference data (probes strategy).

print(predictions_prob)
#>     samp1     samp2     samp3     samp4     samp5 
#> 0.5738002 0.6287368 0.6053712 0.6392831 0.5077835
```

### Count Surrogates

For count outcomes (e.g., cell counts):

``` r
wts_count <- setNames(wts_df$wt_cnt, rownames(wts_df))

predictions_count <- beta_matrix_comp |>
  surro_set(weights = wts_count, intercept = "Intercept") |>
  reference_fill(reference = ref_vec, type = "probes") |>
  surro_calc(transform = "count")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.
#> Filled 25 values using reference data (probes strategy).

print(predictions_count)
#>    samp1    samp2    samp3    samp4    samp5 
#> 2.268177 2.469781 2.485737 2.467513 2.387315
```

## Working with M-values

Some analyses use M-values instead of beta values. The package provides
conversion functions:

``` r
# Load M-value data
data(mval_matrix_comp)

# Convert M-values to beta values
beta_from_m <- convert_m_to_beta(mval_matrix_comp)

# Or convert beta values to M-values
m_from_beta <- convert_beta_to_m(beta_matrix_comp)

# Verify round-trip conversion
all.equal(beta_matrix_comp, 
          convert_m_to_beta(convert_beta_to_m(beta_matrix_comp)))
#> [1] TRUE
```

## Advanced Features

### Getting Diagnostic Information

``` r
# Get detailed diagnostics
results <- surro_calc(
  my_surro,
  transform = "linear",
  return_diagnostics = TRUE
)

# View diagnostics
str(results$diagnostics)
#> List of 14
#>  $ transform                  : chr "linear"
#>  $ n_samples_original         : int 5
#>  $ n_probes_original          : int 10
#>  $ n_weights                  : int 10
#>  $ intercept_used             : logi TRUE
#>  $ n_common_probes            : int 10
#>  $ n_missing_in_weights       : int 0
#>  $ n_missing_in_methyl        : int 0
#>  $ n_completely_missing_probes: int 0
#>  $ n_samples_with_missing     : int 0
#>  $ probes_set_to_zero         : int 0
#>  $ samples_removed            : int 0
#>  $ n_final_samples            : int 5
#>  $ final_result_range         : Named num [1:2] 1.2 1.21
#>   ..- attr(*, "names")= chr [1:2] "min" "max"
```

### Getting Imputation Statistics

``` r
# Get detailed imputation statistics
surro_with_stats <- impute_obs(
  surro_miss,
  method = "mean",
  return_stats = TRUE
)
#> Imputed 4 values in 2 probes using mean method.

# View statistics
print(surro_with_stats$imputation_stats)
#> Imputation Statistics
#> =====================
#> Method: mean
#> Threshold: 0.0% non-missing data required
#> Date: 2025-11-10 20:46:13.503764
#> 
#> Probe Summary:
#> - Total probes: 10
#> - Complete probes: 3 (30.0%)
#> - Partially missing: 2 (20.0%)
#> - Completely missing: 5 (50.0%)
#> 
#> Imputation Results:
#> - Probes imputed: 2 of 2 eligible (100.0%)
#> - Values imputed: 4 of 29 missing (13.8%)
#> - Remaining missing: 25 values
```

### Verbose Output

For detailed information during processing:

``` r
# Verbose mode for debugging
predictions_verbose <- surro_calc(
  my_surro,
  transform = "linear",
  verbose = TRUE
)
```

## Best Practices

1.  **Always check for missing data** using
    [`methyl_miss()`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
    before calculation

2.  **Use appropriate reference values** when filling missing probes
    (ideally from the same platform and tissue type)

3.  **Set reasonable thresholds** for imputation (e.g.,
    `min_nonmiss_prop = 0.5`)

4.  **Choose the correct transformation** based on how the surrogate was
    developed

5.  **Validate results** by comparing to known values or expected ranges

6.  **Document your workflow** including software versions and
    parameters

## Common Issues and Solutions

### Issue: “No common probes found”

**Solution:** Ensure your CpG names match between methylation data and
weights. Check for:

- Different naming conventions (e.g., “cg12345678” vs “cg12345678_1”)
- Leading/trailing whitespace
- Case sensitivity

### Issue: Many samples removed due to missing data

**Solution:**

- Use
  [`impute_obs()`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
  before
  [`surro_calc()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md)
- Lower `min_nonmiss_prop` threshold if appropriate
- Use reference values to fill missing probes

### Issue: Predictions outside expected range

**Solution:**

- Verify you’re using the correct transformation
- Check that input data is on the correct scale (beta vs M-values)
- Ensure reference values are appropriate for your data

## Summary

MethylSurroGetR provides a complete workflow for:

- Organizing methylation data and surrogate weights
- Handling missing probes and observations
- Calculating predictions with various transformations
- Diagnosing and documenting the analysis process

For more information, see the function documentation
([`?surro_calc`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md),
etc.) or visit the [package
website](https://github.com/jagoode27/MethylSurroGetR).

## Session Information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] MethylSurroGetR_0.0.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.37     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.54         cachem_1.1.0      knitr_1.50        htmltools_0.5.8.1
#>  [9] rmarkdown_2.30    lifecycle_1.0.4   cli_3.6.5         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
#> [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        evaluate_1.0.5   
#> [21] bslib_0.9.0       yaml_2.3.10       jsonlite_2.0.0    rlang_1.1.6      
#> [25] fs_1.6.6
```
