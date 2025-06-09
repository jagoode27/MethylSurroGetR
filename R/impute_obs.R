#' Impute Missing Observations in methyl_surro Object
#'
#' This function imputes missing values in the methylation matrix of a \code{methyl_surro} object using either mean or median imputation.
#' It applies row-wise imputation based on a specified non-missing data proportion threshold.
#'
#' @param methyl_surro An object of class \code{methyl_surro}, containing a methylation matrix.
#' @param method A character string indicating the imputation method. Must be either \code{"mean"} or \code{"median"}.
#' @param min_nonmiss_prop The minimum proportion of non-missing data required in a probe (row) for the imputation to proceed. Must be a numeric value between 0 and 1.
#'
#' @return A \code{methyl_surro} object with missing observations imputed, updating the methylation matrix as per the given method and threshold.
#' @importFrom stats median
#'
#' @examples
#' # Load the sample data
#' data(methyl_surro_miss)
#'
#' # Apply mean imputation with a specified threshold for non-missing data in probes
#' result <- impute_obs(methyl_surro_miss, "mean", min_nonmiss_prop = 0.5)
#'
#' # Check the imputed result
#' print(result$methyl)
#'
#' @export
impute_obs <- function(methyl_surro,
                       method = c("mean", "median"),
                       min_nonmiss_prop = 0) {
  method <- match.arg(method)

  if (!is.numeric(min_nonmiss_prop) || min_nonmiss_prop < 0 || min_nonmiss_prop > 1) {
    stop("min_nonmiss_prop must be a numeric value between 0 and 1.")
  }

  methyl <- methyl_surro$methyl
  completely_missing_probes <- rowSums(is.na(methyl)) == ncol(methyl)
  not_imputed_count <- 0

  for (i in seq_len(nrow(methyl))) {
    if (!completely_missing_probes[i]) {
      probe_nonmissing_proportion <- sum(!is.na(methyl[i, ])) / ncol(methyl)

      if (probe_nonmissing_proportion >= min_nonmiss_prop) {
        for (j in seq_len(ncol(methyl))) {
          if (is.na(methyl[i, j])) {
            if (method == "mean") {
              row_mean <- mean(methyl[i, ], na.rm = TRUE)
              methyl[i, j] <- row_mean
            } else if (method == "median") {
              row_median <- median(methyl[i, ], na.rm = TRUE)
              methyl[i, j] <- row_median
            }
          }
        }
      } else {
        # Increment the count for rows not meeting the non-missing proportion threshold
        not_imputed_count <- not_imputed_count + 1
      }
    }
  }

  # Message to inform the user about non-imputed probes
  if (not_imputed_count > 0) {
    message(sprintf("%d probes were not imputed because they did not meet the threshold. Please use methyl_miss() to see details.",
                    not_imputed_count))
  }

  methyl_surro$methyl <- methyl
  return(methyl_surro)
}
