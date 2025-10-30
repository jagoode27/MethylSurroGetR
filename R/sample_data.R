#' Sample Methylation Data for MethylSurroGetR Package
#'
#' This file contains sample datasets used in examples and testing for the MethylSurroGetR package.
#' All datasets are designed to demonstrate package functionality with realistic but manageable data sizes.
#'
#' @section Methylation Data Matrices:
#' The package includes four methylation data matrices demonstrating different data scenarios:
#' \describe{
#'   \item{\code{\link{beta_matrix_comp}}}{Complete beta values without missing data}
#'   \item{\code{\link{beta_matrix_miss}}}{Beta values with strategic missing data patterns}
#'   \item{\code{\link{mval_matrix_comp}}}{Complete M-values derived from beta_matrix_comp}
#'   \item{\code{\link{mval_matrix_miss}}}{M-values with missing data derived from beta_matrix_miss}
#' }
#'
#' @name sample_data
#' @docType data
#' @keywords datasets
NULL

#' Complete Beta Matrix without Missing Values
#'
#' A 15 x 5 matrix of methylation beta values representing complete data without missing observations.
#'
#' @format A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):
#' \describe{
#'   \item{Rows}{CpG probe identifiers (e.g., "cg01", "cg02", etc.)}
#'   \item{Columns}{Sample identifiers ("samp1" through "samp5")}
#' }
#'
#' @source Simulated data
#' @seealso \code{\link{beta_matrix_miss}}, \code{\link{mval_matrix_comp}}, \code{\link{convert_beta_to_m}}
#'
#' @examples
#' data(beta_matrix_comp)
#' str(beta_matrix_comp)
#'
#' @keywords datasets
"beta_matrix_comp"

#' Beta Matrix with Missing Values
#'
#' A 15 x 5 matrix of methylation beta values containing missing values.
#'
#' @format A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):
#' \describe{
#'   \item{Rows}{CpG probe identifiers matching \code{\link{beta_matrix_comp}}}
#'   \item{Columns}{Sample identifiers matching \code{\link{beta_matrix_comp}}}
#' }
#'
#' @source Derived from \code{\link{beta_matrix_comp}} with simulated missing data patterns
#' @seealso \code{\link{beta_matrix_comp}}, \code{\link{impute_obs}}, \code{\link{methyl_miss}}, \code{\link{reference_fill}}
#'
#' @examples
#' data(beta_matrix_miss)
#' str(beta_matrix_miss)
#'
#' @keywords datasets
"beta_matrix_miss"

#' Complete M-Values Matrix
#'
#' A 15 x 5 matrix of methylation M-values representing complete data without missing observations.
#'
#' @format A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):
#' \describe{
#'   \item{Rows}{CpG probe identifiers matching \code{\link{beta_matrix_comp}}}
#'   \item{Columns}{Sample identifiers matching \code{\link{beta_matrix_comp}}}
#' }
#'
#' @source Converted from \code{\link{beta_matrix_comp}} using \code{convert_beta_to_m()}
#' @seealso \code{\link{beta_matrix_comp}}, \code{\link{mval_matrix_miss}}, \code{\link{convert_beta_to_m}}, \code{\link{convert_m_to_beta}}
#'
#' @examples
#' data(mval_matrix_comp)
#' str(mval_matrix_comp)
#'
#' @keywords datasets
"mval_matrix_comp"

#' M-Values Matrix with Missing Values
#'
#' A 15 x 5 matrix of methylation M-values with the same missing values.
#'
#' @format A numeric matrix with 15 rows (CpG probes) and 5 columns (samples):
#' \describe{
#'   \item{Rows}{CpG probe identifiers matching other methylation matrices}
#'   \item{Columns}{Sample identifiers matching other methylation matrices}
#' }
#'
#' @source Converted from \code{\link{beta_matrix_miss}} using \code{convert_beta_to_m()}
#' @seealso \code{\link{mval_matrix_comp}}, \code{\link{beta_matrix_miss}}, \code{\link{convert_beta_to_m}}
#'
#' @examples
#' data(mval_matrix_miss)
#' str(mval_matrix_miss)
#'
#' @keywords datasets
"mval_matrix_miss"

#' Regression Weights Data Frame
#'
#' A data frame containing regression weights for three different simulated statistical models.
#'
#' @format A data frame with 11 rows and 3 columns:
#' \describe{
#'   \item{wt_lin}{Linear regression weights for continuous outcomes}
#'   \item{wt_prb}{Logistic regression weights for binary outcomes (log-odds scale)}
#'   \item{wt_cnt}{Poisson regression weights for count outcomes (log scale)}
#' }
#'
#' Row names represent 10 CpG probe identifiers plus "Intercept" term.
#'
#' Model types demonstrated:
#' \itemize{
#'   \item Linear: Continuous phenotype prediction
#'   \item Logistic: Binary classification (disease status, etc.)
#'   \item Poisson: Count-based outcomes (cell counts, etc.)
#' }
#'
#' @source Simulated data
#'
#' @examples
#' data(wts_df)
#' str(wts_df)
#'
#' # Extract weights vectors
#' linear_weights <- setNames(wts_df$wt_lin, rownames(wts_df))
#' logistic_weights <- setNames(wts_df$wt_prb, rownames(wts_df))
#' poisson_weights <- setNames(wts_df$wt_cnt, rownames(wts_df))
#'
#' @keywords datasets
"wts_df"

#' Reference Values Data Frame
#'
#' A data frame containing population-level reference values for methylation probes,
#' used for imputing completely missing probes when target data lacks specific CpG sites.
#' Values derived from large-scale population methylation studies.
#'
#' @format A data frame with 20 rows (CpG probes) and 2 columns:
#' \describe{
#'   \item{mean}{Population mean beta values for each probe}
#'   \item{median}{Population median beta values for each probe}
#' }
#'
#' Primary applications:
#' \itemize{
#'   \item Filling completely missing probes with \code{\link{reference_fill}}
#'   \item Providing population baselines for comparison
#'   \item Supporting quality control and outlier detection
#' }
#'
#' @source Simulated data
#' @seealso \code{\link{reference_fill}}
#'
#' @examples
#' data(ref_df)
#' str(ref_df)
#'
#' # Extract columns as named vectors
#' ref_mean <- setNames(ref_df$mean, rownames(ref_df))
#' ref_median <- setNames(ref_df$median, rownames(ref_df))
#'
#' @keywords datasets
"ref_df"
