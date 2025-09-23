#' Create a methyl_surro Object
#'
#' This function takes a methylation matrix and a named vector of weights and returns an object
#' of class \code{methyl_surro}. The function aligns the methylation data with the weights,
#' adding missing probes as needed.
#'
#' @param methyl A matrix of methylation data with CpG sites as row names and samples as column names. All values must be numeric.
#' @param weights A named numeric vector with CpG site names, representing the regression weights.
#' @param intercept An optional name of the intercept term within the weights. If provided,
#'        the specified intercept will be separated from the weights and stored separately.
#'
#' @return An object of class \code{methyl_surro} containing the following components:
#' \describe{
#'   \item{methyl}{The methylation matrix, filtered to include only probes present in the weights,
#'   with additional rows of NAs added for weights not present in the methylation data.}
#'   \item{weights}{A named numeric vector of weights (excluding intercept if specified).}
#'   \item{intercept}{The numeric value of the intercept, if provided. \code{NULL} otherwise.}
#' }
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Validates input format and structure
#'   \item Extracts intercept from weights if specified
#'   \item Filters methylation matrix to include only probes present in weights
#'   \item Adds rows of NA values for weight probes not present in methylation data
#'   \item Returns a properly structured \code{methyl_surro} object
#' }
#'
#' Missing probes (those in weights but not in methylation data) should be handled using
#' \code{\link{reference_fill}} before calculating surrogate values.
#'
#' @examples
#' # Load Methylation Beta Matrix
#' data(beta_matrix_comp, package = "MethylSurroGetR")
#' print(beta_matrix_comp)
#'
#' # Load Weights as Named Vector
#' data(wts_vec_lin, package = "MethylSurroGetR")
#' print(wts_vec_lin)
#'
#' # Build the methyl_surro Object
#' surrogate <- surro_set(methyl = beta_matrix_comp,
#'                        weights = wts_vec_lin,
#'                        intercept = "Intercept")
#' print(surrogate)
#'
#' # Example without intercept
#' surrogate_no_int <- surro_set(methyl = beta_matrix_comp,
#'                               weights = wts_vec_lin[-which(names(wts_vec_lin) == "Intercept")])
#' print(surrogate_no_int)
#'
#' @export
surro_set <- function(methyl, weights, intercept = NULL) {

  # Validate methylation matrix
  if (!is.matrix(methyl)) {
    stop("The methylation data must be a matrix.")
  }

  if (!is.numeric(methyl)) {
    stop("All values in the methylation matrix must be numeric.")
  }

  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Methylation matrix cannot be empty.")
  }

  if (is.null(rownames(methyl))) {
    stop("Methylation matrix must have CpG sites as row names.")
  }

  # Validate weights - now restricted to named vectors only
  if (!is.numeric(weights) || !is.vector(weights) || is.null(names(weights))) {
    stop("Weights must be a named numeric vector with CpG site names.")
  }

  if (length(weights) == 0) {
    stop("Weights vector cannot be empty.")
  }

  # Check for duplicate CpG names in weights
  if (any(duplicated(names(weights)))) {
    stop("Duplicate CpG names found in weights.")
  }

  # Check for missing or empty names in weights
  if (any(is.na(names(weights))) || any(names(weights) == "")) {
    stop("All weights must have valid CpG names (no missing or empty names).")
  }

  # Handle intercept extraction
  weights_vector <- weights  # No conversion needed since it's already a vector
  intercept_val <- NULL

  if (!is.null(intercept)) {
    if (!is.character(intercept) || length(intercept) != 1) {
      stop("Intercept must be a single character string.")
    }

    if (intercept %in% names(weights_vector)) {
      intercept_val <- weights_vector[intercept]
      weights_vector <- weights_vector[names(weights_vector) != intercept]

      # Check that we still have weights after removing intercept
      if (length(weights_vector) == 0) {
        stop("No weights remain after removing intercept.")
      }
    } else {
      stop("The specified intercept is not found in the weights.")
    }
  }

  # Filter methylation matrix to include only those probes in weights
  common_cpg <- intersect(rownames(methyl), names(weights_vector))
  methyl_filtered <- methyl[common_cpg, , drop = FALSE]

  # Add rows with NAs for weights not present in methylation data
  missing_cpg <- setdiff(names(weights_vector), rownames(methyl))

  if (length(missing_cpg) > 0) {
    # More efficient approach for large matrices
    missing_matrix <- matrix(NA_real_,
                             nrow = length(missing_cpg),
                             ncol = ncol(methyl),
                             dimnames = list(missing_cpg, colnames(methyl)))
    methyl_filtered <- rbind(methyl_filtered, missing_matrix)
  }

  # Provide informative messages about data alignment
  n_common <- length(common_cpg)
  n_missing <- length(missing_cpg)
  n_weights <- length(weights_vector)

  if (n_missing > 0) {
    message(sprintf("Added %d missing probes with NA values (%.1f%% of weight probes missing from methylation data).",
                    n_missing, 100 * n_missing / n_weights))
  }

  if (n_common < nrow(methyl)) {
    message(sprintf("Filtered methylation matrix from %d to %d probes to match weights.",
                    nrow(methyl), n_common + n_missing))
  }

  # Create methyl_surro object
  methyl_surro <- list(
    methyl = methyl_filtered,
    weights = weights_vector,
    intercept = intercept_val
  )

  class(methyl_surro) <- "methyl_surro"

  return(methyl_surro)
}
