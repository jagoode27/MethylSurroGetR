# imputation_stats Objects

The `imputation_stats` class stores detailed information about a missing
data imputation operation performed by
[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md).
These objects provide comprehensive statistics about what was imputed,
what was skipped, and why.

## Object Structure

An `imputation_stats` object is a list with class `"imputation_stats"`
containing:

- `method`:

  Character string indicating the imputation method used (`"mean"` or
  `"median"`).

- `min_nonmiss_prop`:

  Numeric value (0 to 1) indicating the minimum proportion of
  non-missing data required for a probe to be imputed.

- `timestamp`:

  POSIXct timestamp of when imputation was performed.

- `n_total_probes`:

  Total number of probes in the methylation matrix.

- `n_complete_probes`:

  Number of probes with no missing values.

- `n_completely_missing_probes`:

  Number of probes missing in all samples.

- `n_partially_missing_probes`:

  Number of probes with some (but not all) missing values.

- `n_probes_imputed`:

  Number of probes that had values imputed.

- `n_probes_skipped`:

  Number of partially missing probes that were skipped because they
  didn't meet the `min_nonmiss_prop` threshold.

- `n_values_imputed`:

  Total number of individual missing values that were imputed.

- `n_missing_before_imputation`:

  Total missing values before imputation.

- `n_missing_after_imputation`:

  Total missing values remaining after imputation.

- `imputation_rate`:

  Proportion of originally missing values that were imputed (0 to 1).

- `probes_imputed`:

  Character vector of CpG probe IDs that were imputed.

- `probes_skipped`:

  Character vector of CpG probe IDs that were skipped.

- `values_imputed_per_probe`:

  Named numeric vector showing the number of values imputed for each
  probe (only present when `return_stats = TRUE`).

- `skipped_reasons`:

  Named character vector explaining why each skipped probe was not
  imputed (only present when `return_stats = TRUE`).

## Creating imputation_stats Objects

Objects are created automatically when calling
[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
with `return_stats = TRUE`:

    # Perform imputation with statistics
    result <- impute_obs(
      my_surro,
      method = "mean",
      min_nonmiss_prop = 0.5,
      return_stats = TRUE
    )

    # Access the statistics
    stats <- result$imputation_stats

## Accessing Information

Components can be accessed using standard list notation:

    # View imputation summary
    stats$n_values_imputed
    stats$imputation_rate

    # View which probes were imputed
    stats$probes_imputed

    # View detailed reasons for skipped probes
    stats$skipped_reasons

## Print Method

The package provides a `print` method that displays a formatted summary:

    print(stats)

This shows an organized view of the imputation operation and its
results.

## See also

[`impute_obs`](https://jagoode27.github.io/MethylSurroGetR/reference/impute_obs.md)
for creating objects and performing imputation,
[`methyl_surro`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_surro.md)
for the input object type,
[`reference_fill_stats`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill_stats.md)
for reference filling statistics
