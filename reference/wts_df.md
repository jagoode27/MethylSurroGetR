# Regression Weights Data Frame

A data frame containing regression weights for three different simulated
statistical models.

## Usage

``` r
wts_df
```

## Format

A data frame with 11 rows and 3 columns:

- wt_lin:

  Linear regression weights for continuous outcomes

- wt_prb:

  Logistic regression weights for binary outcomes (log-odds scale)

- wt_cnt:

  Poisson regression weights for count outcomes (log scale)

Row names represent 10 CpG probe identifiers plus "Intercept" term.

Model types demonstrated:

- Linear: Continuous phenotype prediction

- Logistic: Binary classification (disease status, etc.)

- Poisson: Count-based outcomes (cell counts, etc.)

## Source

Simulated data

## Examples

``` r
data(wts_df)
str(wts_df)
#> 'data.frame':    11 obs. of  3 variables:
#>  $ wt_lin: num  -0.00908 -0.00116 0.00598 -0.00756 0.00122 ...
#>  $ wt_prb: num  0.165 -0.405 -0.116 -0.226 0.315 ...
#>  $ wt_cnt: num  0.0509 0.02584 0.04204 -0.09988 -0.00494 ...

# Extract weights vectors
linear_weights <- setNames(wts_df$wt_lin, rownames(wts_df))
logistic_weights <- setNames(wts_df$wt_prb, rownames(wts_df))
poisson_weights <- setNames(wts_df$wt_cnt, rownames(wts_df))
```
