#' Calculate Predicted Values from methyl_surro Object
#'
#' This function calculates predicted values based on the provided methyl_surro object.
#' If any probes (rows) are completely missing, their values are set to zero.
#' Samples with any missing values are removed from the calculation. The function aligns
#' and multiplies the weights with the methylation matrix, adds the intercept if specified,
#' and applies the specified transformation to the results.
#'
#' @param methyl_surro An object of class \code{methyl_surro}.
#' @param transform A character string specifying the transformation to apply. Can be \code{"linear"}, \code{"count"}, or \code{"probability"}.
#' @param verbose Logical. If \code{TRUE}, detailed progress and diagnostic messages are displayed. Default is \code{FALSE}.
#' @param return_diagnostics Logical. If \code{TRUE}, returns additional diagnostic information. Default is \code{FALSE}.
#'
#' @return If \code{return_diagnostics = FALSE}: A named vector of predicted values with names corresponding to the sample names.
#'         If \code{return_diagnostics = TRUE}: A list containing:
#'         \describe{
#'           \item{predictions}{The named vector of predicted values}
#'           \item{diagnostics}{A list of diagnostic information including missing data handling and calculation details}
#'         }
#'
#' @details
#' The function performs validation and processing steps:
#' \itemize{
#'   \item Validates probe-weight alignment between methylation data and weights
#'   \item Sets completely missing probes to zero
#'   \item Removes samples with any missing values
#'   \item Calculates weighted sums and applies transformations
#'   \item Provides detailed feedback about data processing
#' }
#'
#' @examples
#' # Load sample data
#' data(beta_matrix_comp, package = "MethylSurroGetR")
#' data(wts_df, package = "MethylSurroGetR")
#' 
#' # Create weight vector and methyl_surro object
#' wts_vec_lin <- setNames(wts_df$wt_lin, rownames(wts_df))
#' methyl_surro_comp <- surro_set(beta_matrix_comp, wts_vec_lin, "Intercept")
#'
#' # Basic calculation
#' pred_lin <- surro_calc(methyl_surro_comp, transform = "linear")
#' print(pred_lin)
#'
#' # With diagnostic information
#' result_detailed <- surro_calc(methyl_surro_comp,
#'                               transform = "linear",
#'                               return_diagnostics = TRUE,
#'                               verbose = TRUE)
#' print(result_detailed$diagnostics)
#'
#' @export
surro_calc <- function(methyl_surro,
                       transform = c("linear", "count", "probability"),
                       verbose = FALSE,
                       return_diagnostics = FALSE) {

  transform <- match.arg(transform)

  # Input validation
  if (!inherits(methyl_surro, "methyl_surro")) {
    stop("Input must be an object of class 'methyl_surro'.")
  }

  if (is.null(methyl_surro$methyl) || !is.matrix(methyl_surro$methyl)) {
    stop("methyl_surro object must contain a 'methyl' matrix component.")
  }

  if (is.null(methyl_surro$weights) || !is.numeric(methyl_surro$weights)) {
    stop("methyl_surro object must contain a numeric 'weights' component.")
  }

  methyl <- methyl_surro$methyl
  weights <- methyl_surro$weights
  intercept <- methyl_surro$intercept

  # Check for empty data
  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Methylation matrix cannot be empty.")
  }

  if (length(weights) == 0) {
    stop("Weights vector cannot be empty.")
  }

  # Initialize diagnostics tracking
  diagnostics <- list(
    transform = transform,
    n_samples_original = ncol(methyl),
    n_probes_original = nrow(methyl),
    n_weights = length(weights),
    intercept_used = !is.null(intercept)
  )

  if (verbose) {
    message(sprintf("Starting surrogate calculation with %s transformation...", transform))
    message(sprintf("Input: %d probes x %d samples", nrow(methyl), ncol(methyl)))
    message(sprintf("Weights: %d probe weights available", length(weights)))
    if (!is.null(intercept)) {
      message(sprintf("Intercept: %.6f", intercept))
    }
  }

  # Validate probe-weight alignment
  weight_names <- names(weights)
  methyl_names <- rownames(methyl)

  if (is.null(weight_names)) {
    stop("Weights must be a named vector with probe identifiers.")
  }

  if (is.null(methyl_names)) {
    stop("Methylation matrix must have row names (probe identifiers).")
  }

  # Check alignment between weights and methylation data
  common_probes <- intersect(weight_names, methyl_names)
  missing_in_weights <- setdiff(methyl_names, weight_names)
  missing_in_methyl <- setdiff(weight_names, methyl_names)

  if (length(common_probes) == 0) {
    stop("No common probes found between weights and methylation data. Check probe naming.")
  }

  if (length(missing_in_weights) > 0 && verbose) {
    message(sprintf("Note: %d probes in methylation data are not in weights (will be ignored)",
                    length(missing_in_weights)))
  }

  if (length(missing_in_methyl) > 0 && verbose) {
    message(sprintf("Note: %d weight probes are not in methylation data",
                    length(missing_in_methyl)))
  }

  diagnostics$n_common_probes <- length(common_probes)
  diagnostics$n_missing_in_weights <- length(missing_in_weights)
  diagnostics$n_missing_in_methyl <- length(missing_in_methyl)

  # Filter to common probes and align
  methyl <- methyl[common_probes, , drop = FALSE]
  weights <- weights[common_probes]

  # Handle missing data: set completely missing probes to zero
  completely_missing_probes <- rowSums(is.na(methyl)) == ncol(methyl)
  missing_probes <- sum(completely_missing_probes)
  methyl[completely_missing_probes, ] <- 0

  # Remove samples with any missing values
  missing_samples <- colSums(is.na(methyl)) > 0
  omitted_samples <- sum(missing_samples)
  methyl <- methyl[, !missing_samples, drop = FALSE]

  # Store diagnostics about missing data handling
  diagnostics$n_completely_missing_probes <- missing_probes
  diagnostics$n_samples_with_missing <- omitted_samples
  diagnostics$probes_set_to_zero <- missing_probes
  diagnostics$samples_removed <- omitted_samples

  # Display messages about missing data handling
  if (missing_probes > 0 && omitted_samples > 0) {
    message(sprintf("%d probes set to zero, %d samples omitted. Use reference_fill() and impute_obs().",
                    missing_probes, omitted_samples))
  } else if (missing_probes > 0) {
    message(sprintf("%d probes set to zero. Use reference_fill().", missing_probes))
  } else if (omitted_samples > 0) {
    message(sprintf("%d samples omitted. Use impute_obs().", omitted_samples))
  }

  # Sort the weights to match the order of the rows in the methylation matrix
  sorted_weights <- weights[rownames(methyl)]

  # Matrix multiplication to get the predicted values
  if (verbose) message("Calculating weighted sum...")
  pred_values <- as.numeric(t(methyl) %*% sorted_weights)

  # Add the intercept if specified
  if (!is.null(intercept)) {
    pred_values <- pred_values + intercept
    if (verbose) message(sprintf("Added intercept: %.6f", intercept))
  }

  # Apply the specified transformation
  if (verbose) message(sprintf("Applying %s transformation...", transform))

  if (transform == "count") {
    pred_values <- exp(pred_values)
  } else if (transform == "probability") {
    pred_values <- 1 / (1 + exp(-pred_values))
  }

  # Convert to named vector with sample names
  names(pred_values) <- colnames(methyl)

  diagnostics$n_final_samples <- length(pred_values)
  diagnostics$final_result_range <- c(min = min(pred_values), max = max(pred_values))

  if (verbose) {
    message(sprintf("Calculation completed successfully:"))
    message(sprintf("- Final sample count: %d", length(pred_values)))
    message(sprintf("- Result range: %.4f to %.4f", min(pred_values), max(pred_values)))
  }

  # Return results
  if (return_diagnostics) {
    return(list(
      predictions = pred_values,
      diagnostics = diagnostics
    ))
  } else {
    return(pred_values)
  }
}
