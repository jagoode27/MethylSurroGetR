#' Impute Missing Observations in methyl_surro Object
#'
#' This function imputes missing values in the methylation matrix of a \code{methyl_surro} object using either mean or median imputation.
#' It applies row-wise imputation based on a specified non-missing data proportion threshold.
#'
#' @param methyl_surro An object of class \code{methyl_surro}, containing a methylation matrix.
#' @param method A character string indicating the imputation method. Must be either \code{"mean"} or \code{"median"}.
#' @param min_nonmiss_prop The minimum proportion of non-missing data required in a probe (row) for the imputation to proceed. Must be a numeric value between 0 and 1.
#' @param return_stats Logical. If \code{TRUE}, detailed imputation statistics are added to the returned object. Default is \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, detailed progress messages are displayed. Default is \code{FALSE}.
#'
#' @return A \code{methyl_surro} object with missing observations imputed, updating the methylation matrix as per the given method and threshold.
#'         If \code{return_stats = TRUE}, an additional \code{imputation_stats} component is added containing detailed statistics about the imputation process.
#'
#' @details
#' The function uses vectorized operations for efficient processing of large matrices. Probes are only imputed if they meet the
#' minimum non-missing data threshold and are not completely missing.
#'
#' @importFrom stats median
#'
#' @examples
#' # Load the sample data
#' data(methyl_surro_miss)
#'
#' # Apply mean imputation with a specified threshold
#' result <- impute_obs(methyl_surro_miss, "mean", min_nonmiss_prop = 0.5)
#'
#' # Check the imputed result
#' print(result$methyl)
#'
#' # Get detailed imputation statistics
#' result_with_stats <- impute_obs(methyl_surro_miss, "mean",
#'                                 min_nonmiss_prop = 0.5,
#'                                 return_stats = TRUE,
#'                                 verbose = TRUE)
#' print(result_with_stats$imputation_stats)
#'
#' @export
impute_obs <- function(methyl_surro,
                       method = c("mean", "median"),
                       min_nonmiss_prop = 0,
                       return_stats = FALSE,
                       verbose = FALSE) {

  method <- match.arg(method)

  # Validate inputs
  if (!inherits(methyl_surro, "methyl_surro")) {
    stop("Input must be an object of class 'methyl_surro'.")
  }

  if (is.null(methyl_surro$methyl) || !is.matrix(methyl_surro$methyl)) {
    stop("methyl_surro object must contain a 'methyl' matrix component.")
  }

  if (!is.numeric(min_nonmiss_prop) || min_nonmiss_prop < 0 || min_nonmiss_prop > 1) {
    stop("min_nonmiss_prop must be a numeric value between 0 and 1.")
  }

  # Create copy to avoid modifying original
  if (verbose) message("Creating copy of methyl_surro object")
  result_surro <- list(
    methyl = methyl_surro$methyl,
    weights = methyl_surro$weights,
    intercept = methyl_surro$intercept
  )
  class(result_surro) <- "methyl_surro"

  methyl <- result_surro$methyl

  # Check for empty matrix
  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Methylation matrix cannot be empty.")
  }

  # Check if there's any missing data to impute
  if (!any(is.na(methyl))) {
    if (verbose) message("No missing values found. Returning original object unchanged.")
    return(result_surro)
  }

  # Store original for statistics
  original_missing_count <- sum(is.na(methyl))

  if (verbose) {
    message(sprintf("Starting imputation analysis..."))
    message(sprintf("Found %d missing values across %d probes (%.1f%% of total values missing)",
                    original_missing_count, nrow(methyl),
                    100 * original_missing_count / (nrow(methyl) * ncol(methyl))))
    message(sprintf("Using %s imputation with minimum %.1f%% non-missing data threshold",
                    method, 100 * min_nonmiss_prop))
  }

  # Vectorized approach for better performance
  missing_mask <- is.na(methyl)
  probe_missing_counts <- rowSums(missing_mask)
  probe_nonmissing_props <- (ncol(methyl) - probe_missing_counts) / ncol(methyl)

  # Categorize probes
  complete_probes <- probe_missing_counts == 0
  completely_missing_probes <- probe_missing_counts == ncol(methyl)
  partially_missing_probes <- probe_missing_counts > 0 & probe_missing_counts < ncol(methyl)

  # Identify probes eligible for imputation
  eligible_for_imputation <- partially_missing_probes & (probe_nonmissing_props >= min_nonmiss_prop)

  n_complete <- sum(complete_probes)
  n_completely_missing <- sum(completely_missing_probes)
  n_partially_missing <- sum(partially_missing_probes)
  n_eligible <- sum(eligible_for_imputation)
  n_skipped <- sum(partially_missing_probes & !eligible_for_imputation)

  if (verbose) {
    message(sprintf("Probe analysis:"))
    message(sprintf("- %d probes are complete (%.1f%%)", n_complete, 100 * n_complete / nrow(methyl)))
    message(sprintf("- %d probes are completely missing (%.1f%%)", n_completely_missing, 100 * n_completely_missing / nrow(methyl)))
    message(sprintf("- %d probes have partial missing data (%.1f%%)", n_partially_missing, 100 * n_partially_missing / nrow(methyl)))
    message(sprintf("- %d of %d partial probes meet the %.1f%% threshold for imputation",
                    n_eligible, n_partially_missing, 100 * min_nonmiss_prop))
    if (n_skipped > 0) {
      message(sprintf("- %d probes do not meet threshold", n_skipped))
    }
  }

  # Early return if no probes can be imputed
  if (n_eligible == 0) {
    if (verbose) message("No probes meet the imputation threshold.")
    if (return_stats) {
      result_surro$imputation_stats <- create_imputation_stats(
        method, min_nonmiss_prop, methyl, missing_mask,
        n_complete, n_completely_missing, n_partially_missing,
        0, n_skipped, character(0), rownames(methyl)[partially_missing_probes & !eligible_for_imputation],
        0, original_missing_count, original_missing_count
      )
    }
    return(result_surro)
  }

  # Calculate imputation values efficiently
  if (method == "mean") {
    if (verbose) message("Calculating row-wise means for imputation...")
    imputation_values <- rowMeans(methyl, na.rm = TRUE)
  } else if (method == "median") {
    if (verbose) message("Calculating row-wise medians for imputation...")
    imputation_values <- apply(methyl, 1, median, na.rm = TRUE)
  }

  # Vectorized imputation
  if (verbose) message("Performing imputation...")

  imputed_probes <- character(0)
  values_imputed_per_probe <- numeric(0)
  total_values_imputed <- 0

  for (i in which(eligible_for_imputation)) {
    missing_positions <- missing_mask[i, ]
    n_imputed_this_probe <- sum(missing_positions)

    if (n_imputed_this_probe > 0) {
      methyl[i, missing_positions] <- imputation_values[i]
      imputed_probes <- c(imputed_probes, rownames(methyl)[i])
      values_imputed_per_probe <- c(values_imputed_per_probe, n_imputed_this_probe)
      names(values_imputed_per_probe)[length(values_imputed_per_probe)] <- rownames(methyl)[i]
      total_values_imputed <- total_values_imputed + n_imputed_this_probe
    }
  }

  # Final missing count
  final_missing_count <- sum(is.na(methyl))

  # User feedback
  if (verbose && n_eligible > 0) {
    message(sprintf("Imputation completed successfully:"))
    message(sprintf("- Imputed %d missing values across %d probes",
                    total_values_imputed, length(imputed_probes)))
    if (length(imputed_probes) > 0) {
      avg_per_probe <- mean(values_imputed_per_probe)
      message(sprintf("- Average imputation per probe: %.1f values", avg_per_probe))
    }
    if (n_skipped > 0) {
      message(sprintf("- %d probes (%.1f%%) were not imputed due to insufficient data",
                      n_skipped, 100 * n_skipped / nrow(methyl)))
    }
    if (final_missing_count > 0) {
      message(sprintf("- Remaining missing values: %d (%.1f%% of total)",
                      final_missing_count, 100 * final_missing_count / (nrow(methyl) * ncol(methyl))))
    }
  } else if (!verbose) {
    # Minimal feedback for non-verbose mode
    if (total_values_imputed > 0) {
      message(sprintf("Imputed %d values in %d probes using %s method.",
                      total_values_imputed, length(imputed_probes), method))
    }

    if (n_skipped > 0) {
      message(sprintf("%d probes were not imputed because they did not meet the threshold. Use methyl_miss() to see missing probes.",
                      n_skipped))
    }
  }

  # Recommendations for remaining issues
  if (verbose) {
    if (n_skipped > 0) {
      avg_completeness <- mean(probe_nonmissing_props[partially_missing_probes & !eligible_for_imputation])
      message(sprintf("Recommendation: %d probes were skipped (average %.1f%% completeness). Consider:", n_skipped, 100 * avg_completeness))
      message("- Using reference_fill() to handle missing probes")
      message(sprintf("- Lowering min_nonmiss_prop threshold to %.1f", max(0.1, avg_completeness - 0.1)))
      message("- Checking data quality with methyl_miss()")
    }

    if (n_completely_missing > 0) {
      message(sprintf("Note: %d completely missing probes detected. Use reference_fill() to address these before calculating surrogate values.", n_completely_missing))
    }
  }

  # Update the result object
  result_surro$methyl <- methyl

  # Add imputation statistics if requested
  if (return_stats) {
    skipped_probes <- rownames(methyl)[partially_missing_probes & !eligible_for_imputation]
    skipped_reasons <- paste0(round(100 * probe_nonmissing_props[partially_missing_probes & !eligible_for_imputation], 1),
                              "% < ", round(100 * min_nonmiss_prop, 1), "%")
    names(skipped_reasons) <- skipped_probes

    result_surro$imputation_stats <- create_imputation_stats(
      method, min_nonmiss_prop, methyl, missing_mask,
      n_complete, n_completely_missing, n_partially_missing,
      length(imputed_probes), n_skipped, imputed_probes, skipped_probes,
      total_values_imputed, original_missing_count, final_missing_count,
      values_imputed_per_probe, skipped_reasons
    )
  }

  return(result_surro)
}

# Helper function to create imputation statistics
create_imputation_stats <- function(method, min_nonmiss_prop, methyl, missing_mask,
                                    n_complete, n_completely_missing, n_partially_missing,
                                    n_probes_imputed, n_probes_skipped, probes_imputed, probes_skipped,
                                    n_values_imputed, original_missing, final_missing,
                                    values_imputed_per_probe = NULL, skipped_reasons = NULL) {

  stats <- list(
    method = method,
    min_nonmiss_prop = min_nonmiss_prop,
    timestamp = Sys.time(),
    n_total_probes = nrow(methyl),
    n_complete_probes = n_complete,
    n_completely_missing_probes = n_completely_missing,
    n_partially_missing_probes = n_partially_missing,
    n_probes_imputed = n_probes_imputed,
    n_probes_skipped = n_probes_skipped,
    n_values_imputed = n_values_imputed,
    n_missing_before_imputation = original_missing,
    n_missing_after_imputation = final_missing,
    imputation_rate = if(original_missing > 0) n_values_imputed / original_missing else 0,
    probes_imputed = probes_imputed,
    probes_skipped = probes_skipped
  )

  if (!is.null(values_imputed_per_probe)) {
    stats$values_imputed_per_probe <- values_imputed_per_probe
  }

  if (!is.null(skipped_reasons)) {
    stats$skipped_reasons <- skipped_reasons
  }

  class(stats) <- "imputation_stats"
  return(stats)
}

#' Print Method for imputation_stats Objects
#'
#' @param x An object of class \code{imputation_stats}.
#' @param ... Additional arguments passed to print methods.
#' @export
print.imputation_stats <- function(x, ...) {
  cat("Imputation Statistics\n")
  cat("=====================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Threshold: %.1f%% non-missing data required\n", 100 * x$min_nonmiss_prop))
  cat(sprintf("Date: %s\n\n", x$timestamp))

  cat("Probe Summary:\n")
  cat(sprintf("- Total probes: %d\n", x$n_total_probes))
  cat(sprintf("- Complete probes: %d (%.1f%%)\n",
              x$n_complete_probes, 100 * x$n_complete_probes / x$n_total_probes))
  cat(sprintf("- Partially missing: %d (%.1f%%)\n",
              x$n_partially_missing_probes, 100 * x$n_partially_missing_probes / x$n_total_probes))
  cat(sprintf("- Completely missing: %d (%.1f%%)\n\n",
              x$n_completely_missing_probes, 100 * x$n_completely_missing_probes / x$n_total_probes))

  cat("Imputation Results:\n")
  if (x$n_partially_missing_probes > 0) {
    cat(sprintf("- Probes imputed: %d of %d eligible (%.1f%%)\n",
                x$n_probes_imputed, x$n_partially_missing_probes,
                100 * x$n_probes_imputed / x$n_partially_missing_probes))
  } else {
    cat("- No partially missing probes found\n")
  }

  if (x$n_missing_before_imputation > 0) {
    cat(sprintf("- Values imputed: %d of %d missing (%.1f%%)\n",
                x$n_values_imputed, x$n_missing_before_imputation, 100 * x$imputation_rate))
  }

  if (x$n_missing_after_imputation > 0) {
    cat(sprintf("- Remaining missing: %d values\n", x$n_missing_after_imputation))
  } else {
    cat("- No missing values remaining\n")
  }

  if (x$n_probes_skipped > 0) {
    cat(sprintf("\nSkipped probes: %d (insufficient data)\n", x$n_probes_skipped))
  }

  invisible(x)
}
