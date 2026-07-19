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

- `n_samples`:

  Total number of samples in the methylation matrix.

- `n_complete_probes`:

  Number of probes with no missing values after filling.

- `n_completely_missing_probes`:

  Number of completely missing probes remaining after filling.

- `n_partially_missing_probes`:

  Number of partially missing probes remaining after filling.

- `n_reference_probes`:

  Number of probes available in the reference vector.

- `n_probes_filled_complete`:

  Number of completely missing probes that were filled.

- `n_probes_filled_partial`:

  Number of partially missing probes that had individual observations
  filled.

- `n_values_filled_complete`:

  Number of values filled in completely missing probes.

- `n_values_filled_partial`:

  Number of values filled in partially missing probes.

- `n_total_values_filled`:

  Total number of values filled (complete + partial).

- `n_missing_before_filling`:

  Total missing values before filling.

- `n_missing_after_filling`:

  Total missing values remaining after filling.

- `fill_rate`:

  Proportion of originally missing values that were filled (0 to 1).

- `probes_filled_complete`:

  Character vector of CpG probe IDs for completely missing probes that
  were filled.

- `probes_filled_partial`:

  Character vector of CpG probe IDs for probes that had individual
  observations filled.

- `reference_range`:

  Named numeric vector (`min`, `max`) giving the range of the reference
  values used.

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
    stats$probes_filled_complete

    # View detailed filling counts
    stats$n_values_filled_complete

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
