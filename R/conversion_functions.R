#' Convert M-Values to Beta Values
#'
#' Converts methylation M-values to beta values. Beta values represent methylation
#' levels as proportions (0-1), while M-values are log2 ratios that provide better
#' statistical properties for analysis.
#'
#' @param methyl A matrix of methylation values with CpG sites as row names and samples as column names. All values must be numeric.
#' @param in_place Logical. If \code{TRUE}, the conversion is performed in-place to save memory.
#'        \strong{WARNING: This will permanently modify the original matrix.} Default is \code{FALSE}.
#'
#' @return A matrix of beta values corresponding to the input M-values.
#'
#' @details
#' \strong{Memory Usage:} For large matrices, setting \code{in_place = TRUE} can reduce memory usage by approximately 50 percent but will permanently modify the original data.
#'
#' \strong{Numerical Stability:} The function automatically handles extreme M-values that would result in beta values outside the expected range (0 to 1).
#'
#' @seealso \code{\link{convert_beta_to_m}} for the inverse conversion
#'
#' @examples
#' # Load Methylation M-Values Matrix
#' data(mval_matrix_comp, package = "MethylSurroGetR")
#' print(mval_matrix_comp)
#'
#' # Convert M-values to Beta Values (creates new object)
#' beta_values <- convert_m_to_beta(mval_matrix_comp)
#' print(beta_values)
#'
#' # Convert in-place for memory efficiency (modifies original!)
#' \dontrun{
#' # Make a copy first if you need to preserve original
#' mval_copy <- mval_matrix_comp
#' beta_inplace <- convert_m_to_beta(mval_copy, in_place = TRUE)
#' # mval_copy is now modified and contains beta values
#' }
#'
#' @export
convert_m_to_beta <- function(methyl, in_place = FALSE) {
  # Check if the input is a matrix
  if (!is.matrix(methyl)) {
    stop("Input must be a matrix.")
  }

  # Check if all values are numeric
  if (!is.numeric(methyl)) {
    stop("All values in the methylation matrix must be numeric.")
  }

  # Check for empty input
  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Input matrix cannot be empty.")
  }

  # Check for row names (CpG identifiers)
  if (is.null(rownames(methyl))) {
    warning("Matrix lacks row names. Consider adding CpG site identifiers as row names.")
  }

  # Warn user about in_place modification
  if (in_place) {
    message("Using in_place = TRUE will permanently modify the original matrix.")
  }

  # Convert each m-value to a beta value, preserving NAs
  if (in_place) {
    # Store original NA positions
    original_na <- is.na(methyl)

    # Modify original matrix in place
    methyl[] <- 2^methyl / (2^methyl + 1)

    # Restore NAs that were in original data
    methyl[original_na] <- NA

    beta_values <- methyl
  } else {
    # Create new object
    beta_values <- 2^methyl / (2^methyl + 1)
  }

  # Check for values outside expected range and warn
  if (any(!is.na(beta_values) & (beta_values < 0 | beta_values > 1))) {
    warning("Some converted beta values are outside the expected [0,1] range. This may indicate invalid M-values in the input.")
  }

  return(beta_values)
}

#' Convert Beta Values to M-Values
#'
#' Converts methylation beta values to M-values. Beta values represent methylation
#' levels as proportions (0-1), while M-values are log2 ratios that provide better
#' statistical properties for analysis.
#'
#' @param methyl A matrix of methylation values with CpG sites as row names and samples as column names. All values must be numeric.
#' @param in_place Logical. If \code{TRUE}, the conversion is performed in-place to save memory.
#'        \strong{WARNING: This will permanently modify the original matrix.} Default is \code{FALSE}.
#'
#' @return A matrix of M-values corresponding to the input beta values.
#'
#' @details
#' \strong{Memory Usage:} For large matrices, setting \code{in_place = TRUE} can reduce memory usage by approximately 50 percent
#' but will permanently modify the original data.
#'
#' \strong{Numerical Stability:} The function automatically handles extreme beta values (0 or 1) that would otherwise
#' result in infinite M-values by clamping them to a small epsilon away from the boundaries.
#'
#' \strong{Value Ranges:} Beta values should be between 0 and 1. Values outside this range will generate a warning
#' but the conversion will proceed.
#'
#' @seealso \code{\link{convert_m_to_beta}} for the inverse conversion
#'
#' @examples
#' # Load Methylation Beta Matrix
#' data(beta_matrix_comp, package = "MethylSurroGetR")
#' print(beta_matrix_comp)
#'
#' # Convert Beta Values to M-values
#' m_values <- convert_beta_to_m(beta_matrix_comp)
#' print(m_values)
#'
#' # Convert in-place for memory efficiency (modifies original!)
#' \dontrun{
#' beta_copy <- beta_matrix_comp
#' m_inplace <- convert_beta_to_m(beta_copy, in_place = TRUE)
#' # beta_copy is now modified and contains M-values
#' }
#'
#' @export
convert_beta_to_m <- function(methyl, in_place = FALSE) {
  # Check if the input is a matrix
  if (!is.matrix(methyl)) {
    stop("Input must be a matrix.")
  }

  # Check if all values are numeric
  if (!is.numeric(methyl)) {
    stop("All values in the methylation matrix must be numeric.")
  }

  # Check for empty input
  if (nrow(methyl) == 0 || ncol(methyl) == 0) {
    stop("Input matrix cannot be empty.")
  }

  # Check for row names (CpG identifiers)
  if (is.null(rownames(methyl))) {
    warning("Matrix lacks row names. Consider adding CpG site identifiers as row names.")
  }

  # Check for values outside expected range and warn
  if (any(!is.na(methyl) & (methyl < 0 | methyl > 1))) {
    warning("Beta values outside [0,1] range detected. Results may be unreliable.")
  }

  # Warn user about in_place modification
  if (in_place) {
    message("Using in_place = TRUE will permanently modify the original matrix.")
  }

  # Apply numerical stability: clamp extreme beta values and convert
  if (in_place) {
    # Store original NA positions
    original_na <- is.na(methyl)

    # First clamp to safe range
    methyl[] <- pmax(pmin(methyl, 1 - 1e-6), 1e-6)

    # Then convert
    methyl[] <- log2(methyl / (1 - methyl))

    # Restore NAs that were in original data
    methyl[original_na] <- NA

    m_values <- methyl
  } else {
    # Create new object
    # Clamp beta values to prevent infinite M-values
    beta_clamped <- pmax(pmin(methyl, 1 - 1e-6), 1e-6)
    # Convert each beta value to an m-value, preserving NAs
    m_values <- log2(beta_clamped / (1 - beta_clamped))

    # Restore NAs that were in original data
    m_values[is.na(methyl)] <- NA
  }

  return(m_values)
}
