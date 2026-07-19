# Complete methyl_surro Object

A `methyl_surro` object built from complete methylation data, with all
weight probes present. Structurally missing probes have been filled from
reference means using
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md),
so the methylation matrix contains no missing values.

## Usage

``` r
methyl_surro_comp
```

## Format

An object of class `methyl_surro`: a list with three components:

- methyl:

  A 10 x 5 numeric matrix (weight probes x samples) with no missing
  values.

- weights:

  A named numeric vector of 10 linear regression weights (intercept
  removed).

- intercept:

  The numeric intercept value.

## Source

Built from
[`beta_matrix_comp`](https://jagoode27.github.io/MethylSurroGetR/reference/beta_matrix_comp.md),
[`wts_df`](https://jagoode27.github.io/MethylSurroGetR/reference/wts_df.md),
and
[`ref_df`](https://jagoode27.github.io/MethylSurroGetR/reference/ref_df.md)
via
[`surro_set()`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md)
and
[`reference_fill()`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)

## See also

[`methyl_surro_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_surro_miss.md),
[`surro_set`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md),
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md),
[`surro_calc`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md)

## Examples

``` r
data(methyl_surro_comp)
str(methyl_surro_comp)
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
