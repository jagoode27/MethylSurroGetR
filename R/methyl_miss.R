#' Summarize Missing Data in methyl_surro Object
#'
#' This function summarizes the missing data in a methyl_surro object.
#' It identifies probes with partial missing data (missing observations) and probes missing
#' from all samples (missing probes), along with comprehensive summary statistics.
#'
#' @param methyl_surro An object of class \code{methyl_surro}.
#'
#' @return A list of class \code{methyl_miss} with the following elements:
#' \describe{
#'   \item{missing_obs}{A named numeric vector with the proportion of missing observations for each probe where this value is greater than 0 but less than 1.}
#'   \item{missing_probes}{A character vector of probes that are missing in all samples.}
#'   \item{summary}{A list containing summary statistics about the missing data patterns.}
#' }
#'
#' The \code{summary} component contains:
#' \describe{
#'   \item{total_probes}{Total number of probes in the methylation matrix.}
#'   \item{total_samples}{Total number of samples in the methylation matrix.}
#'   \item{n_complete_probes}{Number of probes with no missing values.}
#'   \item{n_missing_obs}{Number of probes with partial missing data.}
#'   \item{n_missing_probes}{Number of completely missing probes.}
#'   \item{overall_missing_rate}{Proportion of all matrix values that are missing.}
#'   \item{missing_obs_rate}{Proportion of probes that have partial missing data.}
#'   \item{missing_probes_rate}{Proportion of probes that are completely missing.}
#'   \item{complete_probes_rate}{Proportion of probes that are complete.}
#' }
#'
#' @examples
#' # Load Methylation Beta Matrix
#' data(beta_matrix_miss, package = "MethylSurroGetR")
#'
#' # Load Weights from a Named Vector
#' data(wts_vec_lin, package = "MethylSurroGetR")
#'
#' # Build the methyl_surro Object
#' surrogate <- surro_set(methyl = beta_matrix_miss,
#'                        weights = wts_vec_lin,
#'                        intercept = "Intercept")
#'
#' # Summarizing Missing Values
#' missing_summary <- methyl_miss(methyl_surro = surrogate)
#' print(missing_summary)
#'
#' # Access specific components
#' missing_summary$missing_obs
#' missing_summary$summary$overall_missing_rate
#'
#' @export
methyl_miss <- function(methyl_surro) {

  # Ensure the input is of class "methyl_surro"
  if (!inherits(methyl_surro, "methyl_surro")) {
    stop("Input must be an object of class 'methyl_surro'.")
  }

  # Validate the methyl component exists and is proper format
  if (is.null(methyl_surro$methyl) || !is.matrix(methyl_surro$methyl)) {
    stop("methyl_surro object must contain a 'methyl' matrix component.")
  }

  methyl <- methyl_surro$methyl

  # Check for empty matrix
  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Methylation matrix cannot be empty.")
  }

  # More efficient calculation for large matrices
  missing_counts <- rowSums(is.na(methyl))
  n_samples <- ncol(methyl)
  n_probes <- nrow(methyl)

  # Vectorized approach
  missing_proportions <- missing_counts / n_samples
  missing_probes <- names(missing_proportions[missing_proportions == 1])
  missing_obs <- missing_proportions[missing_proportions > 0 & missing_proportions < 1]

  # Calculate summary statistics
  n_complete_probes <- sum(missing_counts == 0)
  n_missing_obs <- length(missing_obs)
  n_missing_probes <- length(missing_probes)
  overall_missing_rate <- sum(is.na(methyl)) / (n_probes * n_samples)

  # Create summary statistics
  summary_stats <- list(
    total_probes = n_probes,
    total_samples = n_samples,
    n_complete_probes = n_complete_probes,
    n_missing_obs = n_missing_obs,
    n_missing_probes = n_missing_probes,
    overall_missing_rate = overall_missing_rate,
    missing_obs_rate = n_missing_obs / n_probes,
    missing_probes_rate = n_missing_probes / n_probes,
    complete_probes_rate = n_complete_probes / n_probes
  )

  # Create the result list
  result <- list(
    missing_obs = missing_obs,
    missing_probes = missing_probes,
    summary = summary_stats
  )

  # Assign class for custom print method
  class(result) <- "methyl_miss"

  return(result)
}

#' Print Method for methyl_miss Objects
#'
#' @param x An object of class \code{methyl_miss}.
#' @param ... Additional arguments passed to print methods.
#' @export
print.methyl_miss <- function(x, ...) {
  cat("Missing Data Summary for methyl_surro Object\n")
  cat("============================================\n\n")

  s <- x$summary

  cat(sprintf("Total probes: %d\n", s$total_probes))
  cat(sprintf("Total samples: %d\n", s$total_samples))
  cat(sprintf("Complete probes: %d (%.1f%%)\n",
              s$n_complete_probes,
              100 * s$complete_probes_rate))
  cat(sprintf("Probes with missing observations: %d (%.1f%%)\n",
              s$n_missing_obs,
              100 * s$missing_obs_rate))
  cat(sprintf("Completely missing probes: %d (%.1f%%)\n",
              s$n_missing_probes,
              100 * s$missing_probes_rate))
  cat(sprintf("Overall missing rate: %.1f%%\n\n",
              100 * s$overall_missing_rate))

  if (length(x$missing_obs) > 0) {
    cat("Probes with partial missing data:\n")
    print(head(x$missing_obs, 10))
    if (length(x$missing_obs) > 10) {
      cat(sprintf("... and %d more\n", length(x$missing_obs) - 10))
    }
    cat("\n")
  }

  if (length(x$missing_probes) > 0) {
    cat("Completely missing probes:\n")
    cat(paste(head(x$missing_probes, 10), collapse = ", "))
    if (length(x$missing_probes) > 10) {
      cat(sprintf(", ... and %d more", length(x$missing_probes) - 10))
    }
    cat("\n")
  }

  invisible(x)
}
