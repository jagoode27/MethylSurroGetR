# Beta Matrix with Missing Values

A 15 x 5 matrix of methylation beta values containing missing values.

## Usage

``` r
beta_matrix_miss
```

## Format

A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):

- Rows:

  CpG probe identifiers matching
  [`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md)

- Columns:

  Sample identifiers matching
  [`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md)

## Source

Derived from
[`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md)
with simulated missing data patterns

## See also

[`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md),
[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md),
[`methyl_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md),
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)

## Examples

``` r
data(beta_matrix_miss)
str(beta_matrix_miss)
#>  num [1:15, 1:5] 0.103 0.288 0.435 0.985 0.9 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:15] "cg01" "cg02" "cg04" "cg05" ...
#>   ..$ : chr [1:5] "samp1" "samp2" "samp3" "samp4" ...
```
