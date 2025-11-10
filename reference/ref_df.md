# Reference Values Data Frame

A data frame containing population-level reference values for
methylation probes, used for imputing completely missing probes when
target data lacks specific CpG sites. Values derived from large-scale
population methylation studies.

## Usage

``` r
ref_df
```

## Format

A data frame with 20 rows (CpG probes) and 2 columns:

- mean:

  Population mean beta values for each probe

- median:

  Population median beta values for each probe

Primary applications:

- Filling completely missing probes with
  [`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)

- Providing population baselines for comparison

- Supporting quality control and outlier detection

## Source

Simulated data

## See also

[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)

## Examples

``` r
data(ref_df)
str(ref_df)
#> 'data.frame':    20 obs. of  2 variables:
#>  $ mean  : num  0.399 0.662 0.495 0.528 0.677 ...
#>  $ median: num  0.369 0.788 0.528 0.435 0.782 ...

# Extract columns as named vectors
ref_mean <- setNames(ref_df$mean, rownames(ref_df))
ref_median <- setNames(ref_df$median, rownames(ref_df))
```
