# Convert Beta Values to M-Values

Converts methylation beta values to M-values. Beta values represent
methylation levels as proportions (0-1), while M-values are log2 ratios
that provide better statistical properties for analysis.

## Usage

``` r
convert_beta_to_m(methyl, in_place = FALSE)
```

## Arguments

- methyl:

  A matrix of methylation values with CpG sites as row names and samples
  as column names. All values must be numeric.

- in_place:

  Logical. If `TRUE`, the conversion is performed in-place to save
  memory. **WARNING: This will permanently modify the original matrix.**
  Default is `FALSE`.

## Value

A matrix of M-values corresponding to the input beta values.

## Details

**Memory Usage:** For large matrices, setting `in_place = TRUE` can
reduce memory usage by approximately 50\\ but will permanently modify
the original data.

**Numerical Stability:** The function automatically handles extreme beta
values (0 or 1) that would otherwise result in infinite M-values by
clamping them to a small epsilon away from the boundaries.

**Value Ranges:** Beta values should be between 0 and 1. Values outside
this range will generate a warning but the conversion will proceed.

## See also

[`convert_m_to_beta`](https://jagoode27.github.io/MethylSurroGetR/reference/convert_m_to_beta.md)
for the inverse conversion

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

# Convert Beta Values to M-values
m_values <- convert_beta_to_m(beta_matrix_comp)
print(m_values)
#>           samp1       samp2       samp3      samp4      samp5
#> cg01 -3.1245785 -1.08498899 -0.09870499  1.3668399 -0.7709910
#> cg02 -1.3087821  1.89676790 -0.53119503  2.9161480  3.9816231
#> cg04 -0.3778651 -2.11366803  3.02147005 -2.5916049  5.9627377
#> cg05  6.0328946  1.84533228  3.41784669  0.2853374 -2.4554883
#> cg07  3.1671215 -1.61522389 -4.50943124 -1.0352844  4.3909280
#> cg08  3.0095254  1.17328380  0.83324949  7.4389022  0.9294068
#> cg09  3.0618213 -3.27565268  0.63766815  4.3772858 -3.3195746
#> cg10  2.9649848 -0.19199370 -0.52097826  0.4981981 -2.5961894
#> cg12 -2.2365131  0.06640720 -2.53564040 -0.5579002  1.1543760
#> cg13  4.7029201  3.20716109  1.15908804  1.9594722 -5.3084412
#> cg14 -2.7336497  0.58489613  3.85358688  0.8797447  0.7017176
#> cg16  0.9127963 -1.00331091 -1.21395975 -1.0886523  3.0369612
#> cg17 -2.5856356 -0.49802163 -0.50290941 -0.7749769 -2.4750210
#> cg19 -0.9343795 -0.06572304 -3.95129714 -1.1697374  1.0413099
#> cg20  0.9361366  4.38993756  4.18033193 -1.8279249  1.4871805

# Convert in-place for memory efficiency (modifies original!)
if (FALSE) { # \dontrun{
beta_copy <- beta_matrix_comp
m_inplace <- convert_beta_to_m(beta_copy, in_place = TRUE)
# beta_copy is now modified and contains M-values
} # }
```
