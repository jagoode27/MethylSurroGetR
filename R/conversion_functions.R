#' Conversion Functions for M-Values and Beta Values
#'
#' These functions convert methylation M-values to beta values and vice versa.
#'
#' @param methyl A matrix or data frame of methylation values. All values must be numeric and the matrix or data frame must have CpG sites as row names.
#'
#' @return For `convert_m_to_beta`, the function returns a matrix or data frame of beta values corresponding to the input M-values.
#'
#' For `convert_beta_to_m`, the function returns a matrix or data frame of M-values corresponding to the input beta values.
#'
#' @examples
#' # Load Methylation M-Values Matrix
#' data(mval_matrix_comp, package = "MethylSurroGetR")
#' print(mval_matrix_comp)
#'
#' # Convert M-values to Beta Values
#' beta_values <- convert_m_to_beta(mval_matrix_comp)
#' print(beta_values)
#'
#' # Load Methylation Beta Matrix
#' data(beta_matrix_comp, package = "MethylSurroGetR")
#' print(beta_matrix_comp)
#'
#' # Convert Beta Values to M-values
#' m_values <- convert_beta_to_m(beta_matrix_comp)
#' print(m_values)
#'
#' @name methylConversion
NULL

#' Convert M-Values to Beta Values
#'
#' @rdname methylConversion
#' @export
convert_m_to_beta <- function(methyl) {
  # Check if the input is a matrix or data frame
  if (!is.matrix(methyl) && !is.data.frame(methyl)) {
    stop("Input must be a matrix or data frame.")
  }

  # Preserve input format for return
  is_data_frame <- is.data.frame(methyl)

  # Convert to a consistent format for processing
  if (is_data_frame) {
    m_values <- as.matrix(methyl)
  } else {
    m_values <- methyl
  }

  # Convert each m-value to a beta value, preserving NAs
  beta_values <- 2^m_values / (2^m_values + 1)

  # Convert back to data frame if original input was a data frame
  if (is_data_frame) {
    beta_values <- as.data.frame(beta_values)
  }

  return(beta_values)
}

#' Convert Beta Values to M-Values
#'
#' @rdname methylConversion
#' @export
convert_beta_to_m <- function(methyl) {
  # Check if the input is a matrix or data frame
  if (!is.matrix(methyl) && !is.data.frame(methyl)) {
    stop("Input must be a matrix or data frame.")
  }

  # Preserve input format for return
  is_data_frame <- is.data.frame(methyl)

  # Convert to a consistent format for processing
  if (is_data_frame) {
    beta_values <- as.matrix(methyl)
  } else {
    beta_values <- methyl
  }

  # Convert each beta value to an m-value, preserving NAs
  m_values <- log2(beta_values / (1 - beta_values))

  # Convert back to data frame if original input was a data frame
  if (is_data_frame) {
    m_values <- as.data.frame(m_values)
  }

  return(m_values)
}
