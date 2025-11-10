# Complete M-Values Matrix

A 15 x 5 matrix of methylation M-values representing complete data
without missing observations.

## Usage

``` r
mval_matrix_comp
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

Converted from
[`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md)
using
[`convert_beta_to_m()`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_beta_to_m.md)

## See also

[`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md),
[`mval_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_miss.md),
[`convert_beta_to_m`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_beta_to_m.md),
[`convert_m_to_beta`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_m_to_beta.md)

## Examples

``` r
data(mval_matrix_comp)
str(mval_matrix_comp)
#>  num [1:15, 1:5] -3.125 -1.309 -0.378 6.033 3.167 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:15] "cg01" "cg02" "cg04" "cg05" ...
#>   ..$ : chr [1:5] "samp1" "samp2" "samp3" "samp4" ...
```
