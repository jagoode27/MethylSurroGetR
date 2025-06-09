#' Create a methyl_surro Object
#'
#' This function takes a methylation matrix and a set of weights (provided as a named vector,
#' a single-column data frame, or a single-column matrix with row names) and returns an object
#' of class \code{methyl_surro}.
#'
#' @param methyl A matrix of methylation data. All values must be numeric and the matrix
#'        must have CpG sites as row names.
#' @param weights A named numeric vector or a single-column data frame or matrix
#'        with CpGs as row names, representing the weights.
#' @param intercept An optional name of the intercept term within the weights. If provided,
#'        the specified intercept will be separated from the weights and stored separately.
#'
#' @return An object of class \code{methyl_surro} containing the following components:
#' \describe{
#'   \item{methyl}{The input methylation matrix, filtered to include only probes present in the weights,
#'   and with additional empty rows added for weights not present in the methylation data.}
#'   \item{weights}{A named numeric vector of weights.}
#'   \item{intercept}{The numeric value of the intercept, if provided. \code{NULL} otherwise.}
#' }
#'
#' @examples
#' # Load Methylation Beta Matrix
#' data(beta_matrix_comp, package = "MethylSurroGetR")
#' print(beta_matrix_comp)
#'
#' # Specify Weights from a Named Vector
#' data(wts_vec_lin, package = "MethylSurroGetR")
#' print(wts_vec_lin)
#'
#' # Build the methyl_surro Object
#' surrogate_vec <- surro_set(methyl = beta_matrix_comp,
#'                            weights = wts_vec_lin,
#'                            intercept = "Intercept")
#' print(surrogate_vec)
#'
#' # Specify Weights from a Data Frame
#' data(wts_df, package = "MethylSurroGetR")
#' print(wts_df)
#'
#' # Build the methyl_surro Object
#' surrogate_df <- surro_set(methyl = beta_matrix_comp,
#'                        weights = wts_df["wt_lin"],
#'                        intercept = "Intercept")
#' print(surrogate_df)
#'
#' # Specify Weights from a Matrix
#' data(wts_mat, package = "MethylSurroGetR")
#' print(wts_mat)
#'
#' # Build the methyl_surro Object
#' surrogate_mat <- surro_set(methyl = beta_matrix_comp,
#'                            weights = wts_mat[, "wt_lin"],
#'                            intercept = "Intercept")
#' print(surrogate_mat)
#'
#' @export
surro_set <- function(methyl, weights, intercept = NULL) {

  # Handle methylation matrix
  if (!is.matrix(methyl)) {
    stop("The methylation data must be a matrix.")
  }
  if (!all(is.numeric(methyl))) {
    stop("All values in the methylation matrix must be numeric.")
  }

  # Handle weights input
  weights_vector <- NULL
  if (is.data.frame(weights) || is.matrix(weights)) {
    if (ncol(weights) != 1) {
      stop("Data frame or matrix of weights must be a single column.")
    }
    if (is.null(rownames(weights))) {
      stop("Data frame or matrix must have CpGs as row names.")
    }

    weights_vector <- as.numeric(weights[, 1])
    names(weights_vector) <- rownames(weights)

  } else if (is.numeric(weights) && !is.null(names(weights))) {
    weights_vector <- weights
  } else {
    stop("Weights must be a named numeric vector, or a single-column data frame or matrix with CpGs as row names.")
  }

  # Handle intercept
  intercept_val <- NULL
  if (!is.null(intercept)) {
    if (intercept %in% names(weights_vector)) {
      intercept_val <- weights_vector[intercept]
      weights_vector <- weights_vector[names(weights_vector) != intercept]
    } else {
      stop("The specified intercept is not found in the weights.")
    }
  }

  # Filter methylation matrix to include only those probes in weights
  common_cpg <- intersect(rownames(methyl), names(weights_vector))
  methyl_filtered <- methyl[common_cpg, , drop = FALSE]

  # Add rows with NAs for weights not present in methyl
  missing_cpg <- setdiff(names(weights_vector), rownames(methyl))
  if (length(missing_cpg) > 0) {
    missing_matrix <- matrix(NA, nrow = length(missing_cpg), ncol = ncol(methyl),
                             dimnames = list(missing_cpg, colnames(methyl)))
    methyl_filtered <- rbind(methyl_filtered, missing_matrix)
  }

  # Combine into a single object and assign class
  methyl_surro <- list(methyl = methyl_filtered, weights = weights_vector, intercept = intercept_val)
  class(methyl_surro) <- "methyl_surro"

  # Return the combined object
  return(methyl_surro)
}
