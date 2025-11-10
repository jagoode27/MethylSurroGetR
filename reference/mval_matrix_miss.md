# M-Values Matrix with Missing Values

A 15 x 5 matrix of methylation M-values with the same missing values.

## Usage

``` r
mval_matrix_miss
```

## Format

A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):

- Rows:

  CpG probe identifiers matching other methylation matrices

- Columns:

  Sample identifiers matching other methylation matrices

## Source

Converted from
[`beta_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_miss.md)
using
[`convert_beta_to_m()`](https://jagoode27.github.io/MethylSurroGetR/reference/methylConversion.md)

## See also

[`mval_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/mval_matrix_comp.md),
[`beta_matrix_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_miss.md),
[`convert_beta_to_m`](https://jagoode27.github.io/MethylSurroGetR/reference/methylConversion.md)

## Examples

``` r
data(mval_matrix_miss)
str(mval_matrix_miss)
#>  num [1:15, 1:5] -3.125 -1.309 -0.378 6.033 3.167 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:15] "cg01" "cg02" "cg04" "cg05" ...
#>   ..$ : chr [1:5] "samp1" "samp2" "samp3" "samp4" ...
```
