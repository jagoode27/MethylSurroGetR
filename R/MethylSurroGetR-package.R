#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' MethylSurroGetR: DNA Methylation Surrogate Calculations
#'
#' @description
#' MethylSurroGetR provides a comprehensive toolkit for calculating DNA methylation
#' surrogate biomarkers from existing studies. The package handles the complete workflow
#' from data preparation through prediction calculation, with robust handling of missing
#' data and flexible transformation options.
#'
#' @section Main Workflow Functions:
#' The typical workflow involves four main steps:
#' \describe{
#'   \item{\code{\link{surro_set}}}{Create a surrogate object by combining methylation
#'   data with surrogate weights}
#'   \item{\code{\link{reference_fill}}}{Fill missing CpG probes using reference values
#'   (e.g., from population means)}
#'   \item{\code{\link{impute_obs}}}{Impute missing observations within samples using
#'   mean or median imputation}
#'   \item{\code{\link{surro_calc}}}{Calculate surrogate predictions with linear,
#'   logistic (probability), or Poisson (count) transformations}
#' }
#'
#' @section Data Assessment Functions:
#' \describe{
#'   \item{\code{\link{methyl_miss}}}{Comprehensive missing data assessment for
#'   methylation matrices}
#' }
#'
#' @section Data Conversion Functions:
#' \describe{
#'   \item{\code{\link{convert_beta_to_m}}}{Convert beta values (0-1 scale) to
#'   M-values (log-ratio scale)}
#'   \item{\code{\link{convert_m_to_beta}}}{Convert M-values back to beta values}
#' }
#'
#' @section Sample Datasets:
#' The package includes example datasets for learning and testing:
#' \describe{
#'   \item{\code{\link{beta_matrix_comp}}}{Complete beta value matrix (15 probes × 5 samples)}
#'   \item{\code{\link{beta_matrix_miss}}}{Beta matrix with missing values}
#'   \item{\code{\link{mval_matrix_comp}}}{Complete M-value matrix}
#'   \item{\code{\link{mval_matrix_miss}}}{M-value matrix with missing values}
#'   \item{\code{\link{wts_df}}}{Example surrogate weights with three transformation types}
#'   \item{\code{\link{ref_df}}}{Reference values for missing probe imputation}
#' }
#'
#' @section Getting Started:
#' To get started with MethylSurroGetR:
#' \enumerate{
#'   \item Load your methylation data (beta or M-values) as a matrix with CpG probes
#'   as rows and samples as columns
#'   \item Obtain surrogate weights from a published study or your own model
#'   \item Follow the basic workflow:
#'   \preformatted{
#'   # Create surrogate object
#'   my_surro <- surro_set(methyl_matrix, weights_vector, intercept = "Intercept")
#'   
#'   # Handle missing probes
#'   my_surro <- reference_fill(my_surro, reference_values)
#'   
#'   # Impute missing observations (if needed)
#'   my_surro <- impute_obs(my_surro, method = "mean")
#'   
#'   # Calculate predictions
#'   predictions <- surro_calc(my_surro, transform = "linear")
#'   }
#'   \item See \code{vignette("MethylSurroGetR")} for detailed examples
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Flexible missing data handling with multiple strategies
#'   \item Support for linear, logistic (probability), and Poisson (count) transformations
#'   \item Comprehensive input validation and informative error messages
#'   \item Detailed diagnostic reporting for transparency
#'   \item Memory-efficient operations for large-scale datasets
#'   \item Extensive test coverage ensuring reliability
#' }
#'
#' @section Authors:
#' Joshua A. Goode \email{jagoode@@umich.edu} (ORCID: 0000-0003-3290-0284)
#'
#' @section See Also:
#' Useful links:
#' \itemize{
#'   \item Package website: \url{https://jagoode27.github.io/MethylSurroGetR/}
#'   \item GitHub repository: \url{https://github.com/jagoode27/MethylSurroGetR}
#'   \item Report bugs: \url{https://github.com/jagoode27/MethylSurroGetR/issues}
#' }
#'
#' @docType package
#' @name MethylSurroGetR-package
#' @aliases MethylSurroGetR
NULL