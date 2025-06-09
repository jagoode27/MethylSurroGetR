#' Complete Beta Matrix without Missing Values
#'
#' A 15 x 5 matrix of methylation beta values for sampled probes without any missing data.
#' (Probes are in rows and samples are in columns.)
#'
#' @format A matrix with 15 rows and 5 columns.
#' @examples
#' data(beta_matrix_comp)
#' str(beta_matrix_comp)
#' @name beta_matrix_comp
NULL

#' Beta Matrix with Missing Values
#'
#' A 15 x 5 matrix of methylation beta values for sampled probes with some missing data.
#' (Probes are in rows and samples are in columns.)
#'
#'
#' @format A matrix with 15 rows and 5 columns.
#' @examples
#' data(beta_matrix_miss)
#' str(beta_matrix_miss)
#' @name beta_matrix_miss
NULL

#' Complete M-Values Matrix without Missing Values and its Associated Elements
#'
#' A 15 x 5 matrix of methylation M-values without any missing data.
#' (Probes are in rows and samples are in columns.)
#'
#' @format A matrix with 15 rows and 5 columns.
#' @examples
#' data(mval_matrix_comp)
#' str(mval_matrix_comp)
#' @name mval_matrix_comp
NULL

#' M-Values Matrix with Missing Values
#'
#' A 15 x 5 matrix of methylation M-values with some missing data.
#' (Probes are in rows and samples are in columns.)
#'
#' @format A matrix with 15 rows and 5 columns.
#' @examples
#' data(mval_matrix_miss)
#' str(mval_matrix_miss)
#' @name mval_matrix_miss
NULL

#' Weights Data Frame
#'
#' A data frame containing regression weights for linear, logistic, and Poisson regression models for 10 probes.
#'
#' @format A data frame with 10 rows and 3 columns:
#' \describe{
#'   \item{wt_lin}{Weights for linear regression}
#'   \item{wt_prb}{Weights for logistic (probability) regression}
#'   \item{wt_cnt}{Weights for Poisson (count) regression}
#' }
#' @examples
#' data(wts_df)
#' str(wts_df)
#' @name wts_df
NULL

#' Weights Matrix
#'
#' A 10 x 3 matrix of regression weights.
#'
#' @format A matrix with 10 rows and 3 columns.
#' \describe{
#'   \item{wt_lin}{Weights for linear regression}
#'   \item{wt_prb}{Weights for logistic (probability) regression}
#'   \item{wt_cnt}{Weights for Poisson (count) regression}
#' }
#' @examples
#' data(wts_mat)
#' str(wts_mat)
#' @name wts_mat
NULL

#' Linear Regression Weights Vector
#'
#' A named vector of weights for the linear regression model.
#'
#' @format A named vector with 10 elements.
#' @examples
#' data(wts_vec_lin)
#' str(wts_vec_lin)
#' @name wts_vec_lin
NULL

#' Logistic Regression Weights Vector
#'
#' A named vector of weights for the logistic (probability) regression model.
#'
#' @format A named vector with 10 elements.
#' @examples
#' data(wts_vec_prb)
#' str(wts_vec_prb)
#' @name wts_vec_prb
NULL

#' Poisson Regression Weights Vector
#'
#' A named vector of weights for the Poisson (count) regression model.
#'
#' @format A named vector with 10 elements.
#' @examples
#' data(wts_vec_cnt)
#' str(wts_vec_cnt)
#' @name wts_vec_cnt
NULL

#' Reference Beta Values
#'
#' A data frame containing the mean and median beta values for 20 reference probes.
#'
#' @format A data frame with 20 rows and 2 columns:
#' \describe{
#'   \item{mean}{Mean beta value across samples}
#'   \item{median}{Median beta value across samples}
#' }
#' @examples
#' data(reference)
#' str(reference)
#' @name reference
NULL
