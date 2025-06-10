#' Summarize Missing Data in methyl_surro Object
#'
#' This function summarizes the missing data in a methyl_surro object.
#' It identifies probes with partial missing data (missing observations) and probes missing
#' from all samples (missing probes).
#'
#' @param methyl_surro An object of class \code{methyl_surro}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{missing_obs}{A named numeric vector with the proportion of missing observations for each probe where this value is greater than 0.}
#'   \item{missing_probes}{A character vector of probes that are missing in all samples}
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
#' @export
methyl_miss <- function(methyl_surro) {

  # Ensure the input is of class "methyl_surro"
  if (!inherits(methyl_surro, "methyl_surro")) {
    stop("Input must be an object of class 'methyl_surro'.")
  }

  methyl <- methyl_surro$methyl

  # Identify missing observations
  missing_obs <- rowSums(is.na(methyl)) / ncol(methyl)
  missing_probes <- names(missing_obs[missing_obs == 1])
  missing_obs <- missing_obs[missing_obs > 0 & missing_obs < 1]

  # Create the result list
  result <- list(
    missing_obs = missing_obs,
    missing_probes = missing_probes
  )

  return(result)
}
