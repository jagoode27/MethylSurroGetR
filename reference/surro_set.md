# Create a methyl_surro Object

This function takes a methylation matrix and a named vector of weights
and returns an object of class `methyl_surro`. The function aligns the
methylation data with the weights, adding missing probes as needed.

## Usage

``` r
surro_set(methyl, weights, intercept = NULL)
```

## Arguments

- methyl:

  A matrix of methylation data with CpG sites as row names and samples
  as column names. All values must be numeric.

- weights:

  A named numeric vector with CpG site names, representing the
  regression weights.

- intercept:

  An optional name of the intercept term within the weights. If
  provided, the specified intercept will be separated from the weights
  and stored separately.

## Value

An object of class `methyl_surro` containing the following components:

- methyl:

  The methylation matrix, filtered to include only probes present in the
  weights, with additional rows of NAs added for weights not present in
  the methylation data.

- weights:

  A named numeric vector of weights (excluding intercept if specified).

- intercept:

  The numeric value of the intercept, if provided. `NULL` otherwise.

## Details

The function performs the following operations:

- Validates input format and structure

- Extracts intercept from weights if specified

- Filters methylation matrix to include only probes present in weights

- Adds rows of NA values for weight probes not present in methylation
  data

- Returns a properly structured `methyl_surro` object

Missing probes (those in weights but not in methylation data) should be
handled using
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
before calculating surrogate values.

## Examples

``` r
# Load Methylation Beta Matrix
data(beta_matrix_comp, package = "MethylSurroGetR")
print(beta_matrix_comp)
#>          samp1      samp2      samp3     samp4      samp5
#> cg01 0.1028646 0.32037324 0.48290240 0.7205963 0.36948887
#> cg02 0.2875775 0.78830514 0.40897692 0.8830174 0.94046728
#> cg04 0.4348927 0.18769112 0.89035022 0.1422943 0.98421920
#> cg05 0.9849570 0.78229430 0.91443819 0.5492847 0.15420230
#> cg07 0.8998250 0.24608773 0.04205953 0.3279207 0.95450365
#> cg08 0.8895393 0.69280341 0.64050681 0.9942698 0.65570580
#> cg09 0.8930511 0.09359499 0.60873498 0.9540912 0.09104400
#> cg10 0.8864691 0.46677904 0.41068978 0.5854834 0.14190691
#> cg12 0.1750527 0.51150546 0.14709469 0.4045103 0.69000710
#> cg13 0.9630242 0.90229905 0.69070528 0.7954674 0.02461368
#> cg14 0.1306957 0.59998896 0.93529980 0.6478935 0.61925648
#> cg16 0.6531019 0.33282354 0.30122890 0.3198206 0.89139412
#> cg17 0.1428000 0.41454634 0.41372433 0.3688455 0.15244475
#> cg19 0.3435165 0.48861303 0.06072057 0.3077200 0.67299909
#> cg20 0.6567581 0.95447383 0.94772694 0.2197676 0.73707774

# Load Weights Data Frame and create vector
data(wts_df, package = "MethylSurroGetR")
wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
print(wts_vec_lin)
#>         cg02         cg03         cg06         cg07         cg08         cg11 
#> -0.009083377 -0.001155999  0.005978497 -0.007562015  0.001218960 -0.005869372 
#>         cg13         cg15         cg17         cg18    Intercept 
#> -0.007449367  0.005066157  0.007900907 -0.002510744  1.211000000 

# Build the methyl_surro Object
surrogate <- surro_set(methyl = beta_matrix_comp,
                       weights = wts_vec_lin,
                       intercept = "Intercept")
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.
print(surrogate)
#> $methyl
#>          samp1     samp2      samp3     samp4      samp5
#> cg02 0.2875775 0.7883051 0.40897692 0.8830174 0.94046728
#> cg07 0.8998250 0.2460877 0.04205953 0.3279207 0.95450365
#> cg08 0.8895393 0.6928034 0.64050681 0.9942698 0.65570580
#> cg13 0.9630242 0.9022990 0.69070528 0.7954674 0.02461368
#> cg17 0.1428000 0.4145463 0.41372433 0.3688455 0.15244475
#> cg03        NA        NA         NA        NA         NA
#> cg06        NA        NA         NA        NA         NA
#> cg11        NA        NA         NA        NA         NA
#> cg15        NA        NA         NA        NA         NA
#> cg18        NA        NA         NA        NA         NA
#> 
#> $weights
#>         cg02         cg03         cg06         cg07         cg08         cg11 
#> -0.009083377 -0.001155999  0.005978497 -0.007562015  0.001218960 -0.005869372 
#>         cg13         cg15         cg17         cg18 
#> -0.007449367  0.005066157  0.007900907 -0.002510744 
#> 
#> $intercept
#> Intercept 
#>     1.211 
#> 
#> attr(,"class")
#> [1] "methyl_surro"

# Example without intercept
surrogate_no_int <- surro_set(methyl = beta_matrix_comp,
                              weights = wts_vec_lin[-which(names(wts_vec_lin) == "Intercept")])
#> Added 5 missing probes with NA values (50.0% of weight probes missing from methylation data).
#> Filtered methylation matrix from 15 to 10 probes to match weights.
print(surrogate_no_int)
#> $methyl
#>          samp1     samp2      samp3     samp4      samp5
#> cg02 0.2875775 0.7883051 0.40897692 0.8830174 0.94046728
#> cg07 0.8998250 0.2460877 0.04205953 0.3279207 0.95450365
#> cg08 0.8895393 0.6928034 0.64050681 0.9942698 0.65570580
#> cg13 0.9630242 0.9022990 0.69070528 0.7954674 0.02461368
#> cg17 0.1428000 0.4145463 0.41372433 0.3688455 0.15244475
#> cg03        NA        NA         NA        NA         NA
#> cg06        NA        NA         NA        NA         NA
#> cg11        NA        NA         NA        NA         NA
#> cg15        NA        NA         NA        NA         NA
#> cg18        NA        NA         NA        NA         NA
#> 
#> $weights
#>         cg02         cg03         cg06         cg07         cg08         cg11 
#> -0.009083377 -0.001155999  0.005978497 -0.007562015  0.001218960 -0.005869372 
#>         cg13         cg15         cg17         cg18 
#> -0.007449367  0.005066157  0.007900907 -0.002510744 
#> 
#> $intercept
#> NULL
#> 
#> attr(,"class")
#> [1] "methyl_surro"
```
