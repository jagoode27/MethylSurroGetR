# methyl_surro Objects

The `methyl_surro` class is the core data structure in MethylSurroGetR.
It stores methylation data aligned with surrogate weights, ready for
calculation of surrogate predictions. Objects of this class are created
by
[`surro_set`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md)
and modified by
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
and
[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md).

## Object Structure

A `methyl_surro` object is a list with class `"methyl_surro"`
containing:

- `methyl`:

  A numeric matrix with CpG sites as rows (with CpG IDs as row names)
  and samples as columns (with sample IDs as column names). Values
  represent methylation levels, typically as beta values (0-1 scale) or
  M-values (log-ratio scale).

- `weights`:

  A named numeric vector where names are CpG site IDs and values are the
  regression coefficients for each probe in the surrogate model. The
  intercept term is stored separately (see below).

- `intercept`:

  A single numeric value representing the intercept term of the
  surrogate model, or `NULL` if no intercept is included.

## Creating methyl_surro Objects

Objects are created using
[`surro_set`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md):

    # Load data
    data(beta_matrix_comp)
    data(wts_df)

    # Create weights vector
    wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))

    # Create methyl_surro object
    my_surro <- surro_set(
      methyl = beta_matrix_comp,
      weights = wts_vec,
      intercept = "Intercept"
    )

## Accessing Components

Components can be accessed using standard list notation:

    # Access methylation matrix
    methyl_data <- my_surro$methyl

    # Access weights
    weight_values <- my_surro$weights

    # Access intercept
    intercept_value <- my_surro$intercept

## Modifying methyl_surro Objects

`methyl_surro` objects can be modified by:

- [`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md):
  Fills missing probes and/or observations with reference values

- [`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md):
  Imputes missing observations using mean or median values

These functions return modified `methyl_surro` objects, optionally with
additional statistics components (see below).

## Optional Statistics Components

When `reference_fill` or `impute_obs` are called with
`return_stats = TRUE`, the returned object may contain additional
components:

- `reference_fill_stats`:

  An object of class `reference_fill_stats` containing detailed
  information about the reference filling operation (see
  [`reference_fill_stats`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill_stats.md)).

- `imputation_stats`:

  An object of class `imputation_stats` containing detailed information
  about the imputation operation (see
  [`imputation_stats`](https://jagoode27.github.io/MethylSurroGetR/reference/imputation_stats.md)).

## Using methyl_surro Objects

Once properly prepared, `methyl_surro` objects are used with
[`surro_calc`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md)
to calculate surrogate predictions:

    # Calculate predictions
    predictions <- surro_calc(my_surro, transform = "linear")

## See also

[`surro_set`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_set.md)
for creating objects,
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
and
[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
for handling missing data,
[`surro_calc`](https://jagoode27.github.io/MethylSurroGetR/reference/surro_calc.md)
for calculating predictions,
[`methyl_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
for assessing missing data patterns
