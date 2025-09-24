#' Fill Missing Probes in methyl_surro Object Using Reference Data
#'
#' This function fills in missing probes in the methylation matrix of a \code{methyl_surro} object based on external reference data.
#' Users can specify if imputation should consider only fully missing rows, partially missing rows, or the entire dataset.
#'
#' @param methyl_surro An object of class \code{methyl_surro}, containing a methylation matrix.
#' @param reference A named numeric vector or a matrix/data frame with row names, containing reference values for probes to impute missing data.
#'        If a matrix or data frame is provided, specify the column using the \code{col_name} argument.
#' @param col_name Optional; character string specifying the column to use if \code{reference} is a matrix or data frame.
#' @param type A character string specifying the scope for filling missing data:
#'   - \code{"probes"}: Fill in only missing probes (completely missing rows).
#'   - \code{"obs"}: Fill in only missing observations (partially missing values).
#'   - \code{"all"}: Fill in missing probes and observations.
#' @param return_stats Logical. If \code{TRUE}, detailed filling statistics are added to the returned object. Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, detailed progress messages are displayed. Default is \code{FALSE}.
#'
#' @return A \code{methyl_surro} object with its methylation matrix updated by filling in the specified missing probes based on the reference data.
#'         If \code{return_stats = TRUE}, an additional \code{reference_fill_stats} component is added containing detailed statistics about the filling process.
#'
#' @details
#' The function uses vectorized operations for efficient processing of large matrices. Reference values are validated to ensure they
#' fall within reasonable ranges for methylation data (0-1 for beta values, typically -10 to 10 for M-values).
#'
#' \strong{Filling Strategy:} The function applies reference values consistently across all samples for each probe, which is
#' appropriate for population-level reference data but may not capture individual-level variation.
#'
#' @examples
#' # Load the sample data
#' data(methyl_surro_miss)
#' data(ref_df)
#'
#' # Apply reference filling using a specific column of the reference data
#' result <- reference_fill(methyl_surro = methyl_surro_miss,
#'                          reference = ref_df,
#'                          col_name = "mean",
#'                          type = "probes",
#'                          verbose = TRUE)
#'
#' # Check the result after filling
#' print(result$methyl)
#'
#' # Get detailed filling statistics
#' result_with_stats <- reference_fill(methyl_surro_miss, ref_df,
#'                                      col_name = "mean", type = "all",
#'                                      return_stats = TRUE, verbose = TRUE)
#' print(result_with_stats$reference_fill_stats)
#'
#' @export
reference_fill <- function(methyl_surro, reference, col_name = NULL,
                           type = c("probes", "obs", "all"),
                           return_stats = FALSE, verbose = FALSE) {

  # Validate inputs
  if (!inherits(methyl_surro, "methyl_surro")) {
    stop("Input must be an object of class 'methyl_surro'.")
  }

  if (is.null(methyl_surro$methyl) || !is.matrix(methyl_surro$methyl)) {
    stop("methyl_surro object must contain a 'methyl' matrix component.")
  }

  # Create copy to avoid modifying original
  methyl_surro <- list(
    methyl = methyl_surro$methyl,
    weights = methyl_surro$weights,
    intercept = methyl_surro$intercept
  )
  class(methyl_surro) <- "methyl_surro"

  methyl <- methyl_surro$methyl

  # Check for empty matrix
  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Methylation matrix cannot be empty.")
  }

  type <- match.arg(type)

  # Handle reference data input validation and conversion
  if (is.data.frame(reference) || is.matrix(reference)) {
    if (is.null(col_name)) {
      stop("col_name must be specified when reference is a matrix or data frame.")
    }
    if (!col_name %in% colnames(reference)) {
      stop(sprintf("Column '%s' not found in reference data. Available columns: %s",
                   col_name, paste(colnames(reference), collapse = ", ")))
    }
    if (is.null(rownames(reference))) {
      stop("Reference matrix/data frame must have row names (probe identifiers).")
    }

    reference_vector <- reference[, col_name]
    if (!is.numeric(reference_vector)) {
      stop(sprintf("Reference column '%s' must be numeric.", col_name))
    }

    named_reference <- setNames(reference_vector, rownames(reference))

  } else if (is.numeric(reference) && !is.null(names(reference))) {
    named_reference <- reference
  } else {
    stop("Reference must be a named numeric vector or a matrix/data frame with row names.")
  }

  # Remove NA values from reference and warn
  na_ref_count <- sum(is.na(named_reference))
  if (na_ref_count > 0) {
    named_reference <- named_reference[!is.na(named_reference)]
    if (verbose) {
      message(sprintf("Removed %d reference probes with NA values.", na_ref_count))
    }
  }

  # Check if reference is empty after removing NAs
  if (length(named_reference) == 0) {
    stop("No valid reference values available after removing NAs.")
  }

  # Validate reference values are in reasonable range for methylation data
  ref_min <- min(named_reference)
  ref_max <- max(named_reference)

  if (ref_min < -15 || ref_max > 15) {
    warning("Reference values outside typical methylation range detected. Please verify data type (beta vs M-values).")
  }

  if (ref_min >= 0 && ref_max <= 1) {
    if (verbose) message("Reference values appear to be beta values (0-1 range).")
  } else if (ref_min >= -10 && ref_max <= 10) {
    if (verbose) message("Reference values appear to be M-values.")
  }

  # Pre-compute missing patterns for efficiency
  missing_mask <- is.na(methyl)
  probe_missing_counts <- rowSums(missing_mask)
  n_samples <- ncol(methyl)

  # Categorize probes
  completely_missing_probes <- probe_missing_counts == n_samples
  partially_missing_probes <- probe_missing_counts > 0 & probe_missing_counts < n_samples
  complete_probes <- probe_missing_counts == 0

  completely_missing_probe_names <- rownames(methyl)[completely_missing_probes]
  partially_missing_probe_names <- rownames(methyl)[partially_missing_probes]

  if (verbose) {
    message(sprintf("Starting reference filling analysis..."))
    message(sprintf("Matrix dimensions: %d probes x %d samples", nrow(methyl), ncol(methyl)))
    message(sprintf("Complete probes: %d (%.1f%%)", sum(complete_probes),
                    100 * sum(complete_probes) / nrow(methyl)))
    message(sprintf("Partially missing probes: %d (%.1f%%)", sum(partially_missing_probes),
                    100 * sum(partially_missing_probes) / nrow(methyl)))
    message(sprintf("Completely missing probes: %d (%.1f%%)", sum(completely_missing_probes),
                    100 * sum(completely_missing_probes) / nrow(methyl)))
    message(sprintf("Reference data available for %d probes", length(named_reference)))
  }

  # Check if there's any missing data to fill
  total_missing <- sum(missing_mask)
  if (total_missing == 0) {
    if (verbose) message("No missing values found. Returning original object unchanged.")
    if (return_stats) {
      methyl_surro$reference_fill_stats <- create_reference_fill_stats(
        type, methyl, named_reference, character(0), character(0), 0, 0, 0, 0
      )
    }
    return(methyl_surro)
  }

  # Initialize tracking variables
  filled_probes <- character(0)
  filled_obs <- character(0)
  values_filled_probes <- 0
  values_filled_obs <- 0

  # Fill based on the type using vectorized operations where possible
  if (type %in% c("probes", "all")) {
    # Fill completely missing probes
    available_for_complete <- intersect(completely_missing_probe_names, names(named_reference))

    if (length(available_for_complete) > 0) {
      if (verbose) {
        message(sprintf("Filling %d completely missing probes with reference values...",
                        length(available_for_complete)))
      }

      for (probe in available_for_complete) {
        methyl[probe, ] <- named_reference[probe]
        filled_probes <- c(filled_probes, probe)
        values_filled_probes <- values_filled_probes + n_samples
      }
    } else if (length(completely_missing_probe_names) > 0 && verbose) {
      message("No reference data available for completely missing probes.")
    }
  }

  if (type %in% c("obs", "all")) {
    # Fill partially missing observations
    available_for_partial <- intersect(partially_missing_probe_names, names(named_reference))

    if (length(available_for_partial) > 0) {
      if (verbose) {
        message(sprintf("Filling missing observations in %d partially missing probes...",
                        length(available_for_partial)))
      }

      for (probe in available_for_partial) {
        missing_indices <- missing_mask[probe, ]
        n_missing_this_probe <- sum(missing_indices)

        if (n_missing_this_probe > 0) {
          methyl[probe, missing_indices] <- named_reference[probe]
          filled_obs <- c(filled_obs, probe)
          values_filled_obs <- values_filled_obs + n_missing_this_probe
        }
      }
    } else if (length(partially_missing_probe_names) > 0 && verbose) {
      message("No reference data available for partially missing probes.")
    }
  }

  # Final missing count
  final_missing_count <- sum(is.na(methyl))
  original_missing_count <- total_missing
  total_values_filled <- values_filled_probes + values_filled_obs

  # User feedback
  if (verbose) {
    if (total_values_filled > 0) {
      message(sprintf("Reference filling completed successfully:"))
      if (length(filled_probes) > 0) {
        message(sprintf("- Filled %d completely missing probes (%d values)",
                        length(filled_probes), values_filled_probes))
      }
      if (length(filled_obs) > 0) {
        message(sprintf("- Filled missing observations in %d probes (%d values)",
                        length(filled_obs), values_filled_obs))
      }
      message(sprintf("- Total values filled: %d (%.1f%% of originally missing)",
                      total_values_filled, 100 * total_values_filled / original_missing_count))
    }

    if (final_missing_count > 0) {
      message(sprintf("- Remaining missing values: %d (%.1f%% of total)",
                      final_missing_count, 100 * final_missing_count / (nrow(methyl) * ncol(methyl))))
    } else {
      message("- No missing values remaining")
    }
  } else if (!verbose && total_values_filled > 0) {
    # Minimal feedback for non-verbose mode
    message(sprintf("Filled %d values using reference data (%s strategy).",
                    total_values_filled, type))
  }

  # Check for still completely missing probes after filling and provide guidance
  still_missing_probes <- rownames(methyl)[rowSums(is.na(methyl)) == ncol(methyl)]
  if (length(still_missing_probes) > 0) {
    num_still_missing <- format(length(still_missing_probes), big.mark = ",")
    if (verbose) {
      unavailable_probes <- setdiff(still_missing_probes, names(named_reference))
      if (length(unavailable_probes) > 0) {
        message(sprintf("Note: %d probes remain completely missing (no reference data available).",
                        length(unavailable_probes)))
        message("Consider using impute_obs() for remaining missing observations or obtaining additional reference data.")
      }
    } else {
      message(sprintf("%s probes remain missing. Use methyl_miss() to see missing probes.",
                      num_still_missing))
    }
  }

  # Update the methyl_surro object
  methyl_surro$methyl <- methyl

  # Add reference filling statistics if requested
  if (return_stats) {
    methyl_surro$reference_fill_stats <- create_reference_fill_stats(
      type, methyl_surro$methyl, named_reference, filled_probes, filled_obs,
      values_filled_probes, values_filled_obs, original_missing_count, final_missing_count
    )
  }

  return(methyl_surro)
}

# Helper function to create reference filling statistics
create_reference_fill_stats <- function(type, filled_methyl, reference_data,
                                        filled_probes, filled_obs,
                                        values_filled_probes, values_filled_obs,
                                        original_missing, final_missing) {

  # Calculate missing patterns from filled data for final statistics
  final_missing_mask <- is.na(filled_methyl)
  probe_missing_counts <- rowSums(final_missing_mask)
  n_samples <- ncol(filled_methyl)

  # Calculate original patterns (reconstruct from filled data and statistics)
  completely_missing <- probe_missing_counts == n_samples
  partially_missing <- probe_missing_counts > 0 & probe_missing_counts < n_samples

  stats <- list(
    type = type,
    timestamp = Sys.time(),
    n_total_probes = nrow(filled_methyl),
    n_samples = n_samples,
    n_complete_probes = sum(probe_missing_counts == 0),
    n_completely_missing_probes = sum(completely_missing),
    n_partially_missing_probes = sum(partially_missing),
    n_reference_probes = length(reference_data),
    n_probes_filled_complete = length(filled_probes),
    n_probes_filled_partial = length(filled_obs),
    n_values_filled_complete = values_filled_probes,
    n_values_filled_partial = values_filled_obs,
    n_total_values_filled = values_filled_probes + values_filled_obs,
    n_missing_before_filling = original_missing,
    n_missing_after_filling = final_missing,
    fill_rate = if(original_missing > 0) (values_filled_probes + values_filled_obs) / original_missing else 0,
    probes_filled_complete = filled_probes,
    probes_filled_partial = filled_obs,
    reference_range = c(min = min(reference_data), max = max(reference_data))
  )

  class(stats) <- "reference_fill_stats"
  return(stats)
}

#' Print Method for reference_fill_stats Objects
#'
#' @param x An object of class \code{reference_fill_stats}.
#' @param ... Additional arguments passed to print methods.
#' @export
print.reference_fill_stats <- function(x, ...) {
  cat("Reference Filling Statistics\n")
  cat("============================\n")
  cat(sprintf("Type: %s\n", x$type))
  cat(sprintf("Date: %s\n\n", x$timestamp))

  cat("Matrix Summary:\n")
  cat(sprintf("- Total probes: %d\n", x$n_total_probes))
  cat(sprintf("- Total samples: %d\n", x$n_samples))
  cat(sprintf("- Complete probes: %d (%.1f%%)\n",
              x$n_complete_probes, 100 * x$n_complete_probes / x$n_total_probes))
  cat(sprintf("- Partially missing: %d (%.1f%%)\n",
              x$n_partially_missing_probes, 100 * x$n_partially_missing_probes / x$n_total_probes))
  cat(sprintf("- Completely missing: %d (%.1f%%)\n\n",
              x$n_completely_missing_probes, 100 * x$n_completely_missing_probes / x$n_total_probes))

  cat("Reference Data:\n")
  cat(sprintf("- Available reference probes: %d\n", x$n_reference_probes))
  cat(sprintf("- Reference value range: %.3f to %.3f\n\n",
              x$reference_range["min"], x$reference_range["max"]))

  cat("Filling Results:\n")
  if (x$n_probes_filled_complete > 0) {
    cat(sprintf("- Complete probes filled: %d (%d values)\n",
                x$n_probes_filled_complete, x$n_values_filled_complete))
  }
  if (x$n_probes_filled_partial > 0) {
    cat(sprintf("- Partial probes filled: %d (%d values)\n",
                x$n_probes_filled_partial, x$n_values_filled_partial))
  }

  if (x$n_missing_before_filling > 0) {
    cat(sprintf("- Total values filled: %d of %d missing (%.1f%%)\n",
                x$n_total_values_filled, x$n_missing_before_filling, 100 * x$fill_rate))
  }

  if (x$n_missing_after_filling > 0) {
    cat(sprintf("- Remaining missing: %d values\n", x$n_missing_after_filling))
  } else {
    cat("- No missing values remaining\n")
  }

  invisible(x)
}
