#' Fill Missing Probes in methyl_surro Object Using Reference Data
#'
#' This function fills in missing probes in the methylation matrix of a \code{methyl_surro} object based on external reference data.
#' Users can specify if imputation should consider only fully missing rows, partially missing rows, or the entire dataset.
#'
#' @param methyl_surro An object of class \code{methyl_surro}, containing a methylation matrix.
#' @param reference A named numeric vector or a matrix/data frame with row names, containing reference values for probes to impute missing data.
#'        If a matrix or data frame is provided, specify the column using the \code{col_name} argument.
#' @param col_name Optional; character string specifying the column to use if \code{reference} is a matrix or data frame.
#' @param type A character string specifying the scope for filling missing data:
#'   - \code{"probes"}: Fill in only missing probes.
#'   - \code{"obs"}: Fill in only missing observations
#'   - \code{"all"}: Fill in missing probes and observations
#'
#' @return A \code{methyl_surro} object with its methylation matrix updated by filling in the specified missing probes based on the reference data.
#'
#' @examples
#' # Load the sample data
#' data(methyl_surro_miss)
#' data(reference)
#'
#' # Apply reference filling using a specific column of the reference data
#' result <- reference_fill(methyl_surro = methyl_surro_miss,
#'                          reference = reference,
#'                          col_name = "mean",
#'                          type = "probes")
#'
#' # Check the result after filling
#' print(result$methyl)
#'
#' @export
reference_fill <- function(methyl_surro, reference, col_name = NULL, type = c("probes", "obs", "all")) {
  # Handle reference data input for data frame or matrix
  if (is.data.frame(reference) || is.matrix(reference)) {
    if (is.null(col_name) || !col_name %in% colnames(reference)) {
      stop("Please provide a valid col_name indicating which column of the reference data frame or matrix to use.")
    }
    reference_vector <- reference[, col_name]
    named_reference <- setNames(reference_vector, rownames(reference))
  } else if (is.numeric(reference) && !is.null(names(reference))) {
    named_reference <- reference
  } else {
    stop("reference must be a named numeric vector or a matrix/data frame with row names, and col_name should be provided if reference is a matrix or data frame.")
  }

  type <- match.arg(type)

  methyl <- methyl_surro$methyl

  # Identify completely missing rows (probes)
  completely_missing_probes <- rownames(methyl)[rowSums(is.na(methyl)) == ncol(methyl)]

  # Fill based on the type
  if (type == "probes") {
    for (probe in completely_missing_probes) {
      if (probe %in% names(named_reference)) {
        methyl[probe, ] <- named_reference[probe]
      }
    }
  } else if (type == "obs") {
    not_completely_missing_probes <- setdiff(rownames(methyl), completely_missing_probes)
    for (probe in not_completely_missing_probes) {
      missing_indices <- is.na(methyl[probe, ])
      if (any(missing_indices) && probe %in% names(named_reference)) {
        methyl[probe, missing_indices] <- named_reference[probe]
      }
    }
  } else if (type == "all") {
    for (probe in rownames(methyl)) {
      missing_indices <- is.na(methyl[probe, ])
      if (any(missing_indices) && probe %in% names(named_reference)) {
        methyl[probe, missing_indices] <- named_reference[probe]
      }
    }
  }

  # Check for still completely missing probes after filling
  still_missing_probes <- rownames(methyl)[rowSums(is.na(methyl)) == ncol(methyl)]
  if (length(still_missing_probes) > 0) {
    num_still_missing <- format(length(still_missing_probes), big.mark = ",")
    message(sprintf("%s probes remain missing. Please use methyl_miss() to see missing probes.",
                    num_still_missing))
  }

  methyl_surro$methyl <- methyl
  return(methyl_surro)
}
