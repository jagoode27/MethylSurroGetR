# reference_fill_stats Objects

The `reference_fill_stats` class stores detailed information about a
reference filling operation performed by
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md).
These objects provide comprehensive statistics about what was filled and
how the missing data was addressed.

## Object Structure

A `reference_fill_stats` object is a list with class
`"reference_fill_stats"` containing:

- `type`:

  Character string indicating the fill type used: `"probes"` (only
  completely missing probes), `"obs"` (only missing observations), or
  `"all"` (both).

- `timestamp`:

  POSIXct timestamp of when filling was performed.

- `n_total_probes`:

  Total number of probes in the methylation matrix.

- `n_total_samples`:

  Total number of samples in the methylation matrix.

- `n_completely_missing_probes_before`:

  Number of completely missing probes before filling.

- `n_completely_missing_probes_after`:

  Number of completely missing probes after filling.

- `n_probes_filled`:

  Number of completely missing probes that were filled.

- `n_obs_filled`:

  Number of individual missing observations that were filled (in
  partially missing probes).

- `n_total_values_filled`:

  Total number of values filled (sum of probes filled × samples +
  individual observations filled).

- `n_missing_before`:

  Total missing values before filling.

- `n_missing_after`:

  Total missing values after filling.

- `fill_rate`:

  Proportion of originally missing values that were filled (0 to 1).

- `probes_filled`:

  Character vector of CpG probe IDs for completely missing probes that
  were filled.

- `obs_filled`:

  Character vector of CpG probe IDs for probes that had individual
  observations filled.

- `values_filled_per_probe`:

  Named numeric vector showing the number of values filled for each
  probe (only present when `return_stats = TRUE`).

- `probes_not_in_reference`:

  Character vector of CpG probe IDs that were missing but not present in
  the reference values (only present when applicable).

## Creating reference_fill_stats Objects

Objects are created automatically when calling
[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
with `return_stats = TRUE`:

    # Perform reference filling with statistics
    result <- reference_fill(
      my_surro,
      reference = ref_values,
      type = "probes",
      return_stats = TRUE
    )

    # Access the statistics
    stats <- result$reference_fill_stats

## Accessing Information

Components can be accessed using standard list notation:

    # View filling summary
    stats$n_total_values_filled
    stats$fill_rate

    # View which probes were filled
    stats$probes_filled

    # View detailed filling counts
    stats$values_filled_per_probe

## Print Method

The package provides a `print` method that displays a formatted summary:

    print(stats)

This shows an organized view of the reference filling operation and its
results.

## See also

[`reference_fill`](https://jagoode27.github.io/MethylSurroGetR/reference/reference_fill.md)
for creating objects and performing reference filling,
[`methyl_surro`](https://jagoode27.github.io/MethylSurroGetR/reference/methyl_surro.md)
for the input object type,
[`imputation_stats`](https://jagoode27.github.io/MethylSurroGetR/reference/imputation_stats.md)
for imputation statistics
