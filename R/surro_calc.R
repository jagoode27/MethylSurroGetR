#' Calculate Predicted Values from methyl_surro Object
#'
#' This function calculates predicted values based on the provided methyl_surro object.
#' If any probes (rows) are completely missing, their values are set to zero, and a count of these
#' probes is stored. Columns (samples) with any missing values are removed, and a count of these
#' is also recorded. The function aligns and multiplies the weights with the methylation matrix,
#' adds the intercept if specified, and applies the specified transformation to the results.
#'
#' @param methyl_surro An object of class \code{methyl_surro}.
#' @param transform A character string specifying the transformation to apply. It can be \code{"linear"}, \code{"count"}, or \code{"probability"}.
#'
#' @return A named vector of predicted values with names corresponding to the sample names in the methylation matrix.
#'
#' @examples
#' # Load \code{methyl_surro} Object
#' data(beta_matrix_comp, package = "MethylSurroGetR")
#'
#' # Generating Surrogate Estimates from Linear Model Weights
#' pred_lin <- surro_calc(methyl_surro = methyl_surro_comp,
#'                        transform = "linear")
#' print(pred_lin)
#'
#' @export
surro_calc <- function(methyl_surro, transform = c("linear", "count", "probability")) {
  transform <- match.arg(transform)

  # Ensure the input is of class "methyl_surro"
  if (!inherits(methyl_surro, "methyl_surro")) {
    stop("Input must be an object of class 'methyl_surro'.")
  }

  methyl <- methyl_surro$methyl
  weights <- methyl_surro$weights
  intercept <- methyl_surro$intercept

  # Replace missing values in completely missing rows with zero
  completely_missing_probes <- rowSums(is.na(methyl)) == ncol(methyl)
  missing_probes <- sum(completely_missing_probes)
  methyl[completely_missing_probes, ] <- 0

  # Remove columns with any missing values
  missing_samples <- colSums(is.na(methyl)) > 0
  omitted_samples <- sum(missing_samples)
  methyl <- methyl[, !missing_samples, drop = FALSE]

  # Display messages
  if (missing_probes > 0 && omitted_samples > 0) {
    message(sprintf("%d probes were set to zero due to missing values, and %d samples were omitted due to missing data.\nIt is strongly recommended that missing probes be filled in with reference_fill() and missing observations be imputed in with impute_obs().",
                    missing_probes, omitted_samples))
  } else if (missing_probes > 0) {
    message(sprintf("%d probes were set to zero due to missing values.\nIt is strongly recommended that missing probes be filled in with reference_fill().", missing_probes))
  } else if (omitted_samples > 0) {
    message(sprintf("%d samples were omitted due to missing data.\nIt is strongly recommended that missing observations be imputed in with impute_obs().", omitted_samples))
  }

  # Sort the weights to match the order of the rows in the methylation matrix
  sorted_weights <- weights[rownames(methyl)]

  # Matrix multiplication to get the predicted values
  pred_values <- as.numeric(t(methyl) %*% sorted_weights)

  # Add the intercept if specified
  if (!is.null(intercept)) {
    pred_values <- pred_values + intercept
  }

  # Apply the specified transformation
  if (transform == "count") {
    pred_values <- exp(pred_values)
  } else if (transform == "probability") {
    pred_values <- 1 / (1 + exp(-pred_values))
  }

  # Convert to named vector with sample names
  names(pred_values) <- colnames(methyl)

  return(pred_values)
}
