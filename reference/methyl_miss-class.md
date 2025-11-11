# methyl_miss Objects

The `methyl_miss` class stores comprehensive information about missing
data patterns in a
[`methyl_surro`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_surro.md)
object. Objects of this class are created by the
[`methyl_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
function and provide detailed summaries of both partial and complete
missingness.

## Object Structure

A `methyl_miss` object is a list with class `"methyl_miss"` containing:

- `missing_obs`:

  A named numeric vector showing the proportion of missing observations
  (0 to 1) for each probe that has partial missingness (i.e., missing in
  some but not all samples). Names are CpG probe IDs.

- `missing_probes`:

  A character vector of CpG probe IDs that are completely missing (i.e.,
  missing in all samples).

- `summary`:

  A list containing overall summary statistics (see details below).

## Summary Statistics

The `summary` component contains:

- `total_probes`:

  Total number of CpG probes in the methylation matrix.

- `total_samples`:

  Total number of samples in the methylation matrix.

- `n_complete_probes`:

  Number of probes with no missing values.

- `n_missing_obs`:

  Number of probes with partial missing data (some but not all samples
  missing).

- `n_missing_probes`:

  Number of completely missing probes (all samples missing).

- `overall_missing_rate`:

  Proportion of all matrix values that are missing (0 to 1).

- `missing_obs_rate`:

  Proportion of probes with partial missing data (0 to 1).

- `missing_probes_rate`:

  Proportion of probes that are completely missing (0 to 1).

- `complete_probes_rate`:

  Proportion of probes with no missing values (0 to 1).

## Creating methyl_miss Objects

Objects are created using
[`methyl_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md):

    # Create methyl_surro object
    my_surro <- surro_set(beta_matrix_miss, weights_vec, "Intercept")

    # Assess missing data
    miss_info <- methyl_miss(my_surro)

## Accessing Information

Components can be accessed using standard list notation:

    # View probes with partial missingness
    miss_info$missing_obs

    # View completely missing probes
    miss_info$missing_probes

    # Access summary statistics
    miss_info$summary$overall_missing_rate
    miss_info$summary$n_missing_probes

## Print Method

The package provides a `print` method that displays a formatted summary:

    print(miss_info)

This shows an organized view of missing data patterns and summary
statistics.

## See also

[`methyl_miss`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_miss.md)
for creating objects,
[`methyl_surro`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_surro.md)
for the input object type,
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
and
[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
for addressing missing data
