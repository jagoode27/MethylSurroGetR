# Complete Beta Matrix without Missing Values

A 15 x 5 matrix of methylation beta values representing complete data
without missing observations.

## Usage

``` r
beta_matrix_comp
```

## Format

A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):

- Rows:

  CpG probe identifiers (e.g., "cg01", "cg02", etc.)

- Columns:

  Sample identifiers ("samp1" through "samp5")

## Source

Simulated data

## See also

[`beta_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_miss.md),
[`mval_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_comp.md),
[`convert_beta_to_m`](https://jagoode27.github.io/MethylSurroGetR/reference/methylConversion.md)

## Examples

``` r
data(beta_matrix_comp)
str(beta_matrix_comp)
#>  num [1:15, 1:5] 0.103 0.288 0.435 0.985 0.9 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:15] "cg01" "cg02" "cg04" "cg05" ...
#>   ..$ : chr [1:5] "samp1" "samp2" "samp3" "samp4" ...
```
