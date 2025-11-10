#' methyl_surro Objects
#'
#' @description
#' The \code{methyl_surro} class is the core data structure in MethylSurroGetR. It stores
#' methylation data aligned with surrogate weights, ready for calculation of surrogate
#' predictions. Objects of this class are created by \code{\link{surro_set}} and modified
#' by \code{\link{reference_fill}} and \code{\link{impute_obs}}.
#'
#' @section Object Structure:
#' A \code{methyl_surro} object is a list with class \code{"methyl_surro"} containing:
#' \describe{
#'   \item{\code{methyl}}{A numeric matrix with CpG sites as rows (with CpG IDs as row names)
#'   and samples as columns (with sample IDs as column names). Values represent methylation
#'   levels, typically as beta values (0-1 scale) or M-values (log-ratio scale).}
#'   \item{\code{weights}}{A named numeric vector where names are CpG site IDs and values
#'   are the regression coefficients for each probe in the surrogate model. The intercept
#'   term is stored separately (see below).}
#'   \item{\code{intercept}}{A single numeric value representing the intercept term of the
#'   surrogate model, or \code{NULL} if no intercept is included.}
#' }
#'
#' @section Creating methyl_surro Objects:
#' Objects are created using \code{\link{surro_set}}:
#' \preformatted{
#' # Load data
#' data(beta_matrix_comp)
#' data(wts_df)
#'
#' # Create weights vector
#' wts_vec <- setNames(wts_df$wt_lin, rownames(wts_df))
#'
#' # Create methyl_surro object
#' my_surro <- surro_set(
#'   methyl = beta_matrix_comp,
#'   weights = wts_vec,
#'   intercept = "Intercept"
#' )
#' }
#'
#' @section Accessing Components:
#' Components can be accessed using standard list notation:
#' \preformatted{
#' # Access methylation matrix
#' methyl_data <- my_surro$methyl
#'
#' # Access weights
#' weight_values <- my_surro$weights
#'
#' # Access intercept
#' intercept_value <- my_surro$intercept
#' }
#'
#' @section Modifying methyl_surro Objects:
#' \code{methyl_surro} objects can be modified by:
#' \itemize{
#'   \item \code{\link{reference_fill}}: Fills missing probes and/or observations with
#'   reference values
#'   \item \code{\link{impute_obs}}: Imputes missing observations using mean or median values
#' }
#'
#' These functions return modified \code{methyl_surro} objects, optionally with additional
#' statistics components (see below).
#'
#' @section Optional Statistics Components:
#' When \code{reference_fill} or \code{impute_obs} are called with \code{return_stats = TRUE},
#' the returned object may contain additional components:
#' \describe{
#'   \item{\code{reference_fill_stats}}{An object of class \code{reference_fill_stats}
#'   containing detailed information about the reference filling operation (see
#'   \code{\link{reference_fill_stats}}).}
#'   \item{\code{imputation_stats}}{An object of class \code{imputation_stats} containing
#'   detailed information about the imputation operation (see \code{\link{imputation_stats}}).}
#' }
#'
#' @section Using methyl_surro Objects:
#' Once properly prepared, \code{methyl_surro} objects are used with \code{\link{surro_calc}}
#' to calculate surrogate predictions:
#' \preformatted{
#' # Calculate predictions
#' predictions <- surro_calc(my_surro, transform = "linear")
#' }
#'
#' @seealso
#' \code{\link{surro_set}} for creating objects,
#' \code{\link{reference_fill}} and \code{\link{impute_obs}} for handling missing data,
#' \code{\link{surro_calc}} for calculating predictions,
#' \code{\link{methyl_miss}} for assessing missing data patterns
#'
#' @name methyl_surro
NULL

#' methyl_miss Objects
#'
#' @description
#' The \code{methyl_miss} class stores comprehensive information about missing data patterns
#' in a \code{\link{methyl_surro}} object. Objects of this class are created by the
#' \code{\link{methyl_miss}} function and provide detailed summaries of both partial and
#' complete missingness.
#'
#' @section Object Structure:
#' A \code{methyl_miss} object is a list with class \code{"methyl_miss"} containing:
#' \describe{
#'   \item{\code{missing_obs}}{A named numeric vector showing the proportion of missing
#'   observations (0 to 1) for each probe that has partial missingness (i.e., missing in
#'   some but not all samples). Names are CpG probe IDs.}
#'   \item{\code{missing_probes}}{A character vector of CpG probe IDs that are completely
#'   missing (i.e., missing in all samples).}
#'   \item{\code{summary}}{A list containing overall summary statistics (see details below).}
#' }
#'
#' @section Summary Statistics:
#' The \code{summary} component contains:
#' \describe{
#'   \item{\code{total_probes}}{Total number of CpG probes in the methylation matrix.}
#'   \item{\code{total_samples}}{Total number of samples in the methylation matrix.}
#'   \item{\code{n_complete_probes}}{Number of probes with no missing values.}
#'   \item{\code{n_missing_obs}}{Number of probes with partial missing data (some but not
#'   all samples missing).}
#'   \item{\code{n_missing_probes}}{Number of completely missing probes (all samples missing).}
#'   \item{\code{overall_missing_rate}}{Proportion of all matrix values that are missing (0 to 1).}
#'   \item{\code{missing_obs_rate}}{Proportion of probes with partial missing data (0 to 1).}
#'   \item{\code{missing_probes_rate}}{Proportion of probes that are completely missing (0 to 1).}
#'   \item{\code{complete_probes_rate}}{Proportion of probes with no missing values (0 to 1).}
#' }
#'
#' @section Creating methyl_miss Objects:
#' Objects are created using \code{\link{methyl_miss}}:
#' \preformatted{
#' # Create methyl_surro object
#' my_surro <- surro_set(beta_matrix_miss, weights_vec, "Intercept")
#'
#' # Assess missing data
#' miss_info <- methyl_miss(my_surro)
#' }
#'
#' @section Accessing Information:
#' Components can be accessed using standard list notation:
#' \preformatted{
#' # View probes with partial missingness
#' miss_info$missing_obs
#'
#' # View completely missing probes
#' miss_info$missing_probes
#'
#' # Access summary statistics
#' miss_info$summary$overall_missing_rate
#' miss_info$summary$n_missing_probes
#' }
#'
#' @section Print Method:
#' The package provides a \code{print} method that displays a formatted summary:
#' \preformatted{
#' print(miss_info)
#' }
#' This shows an organized view of missing data patterns and summary statistics.
#'
#' @seealso
#' \code{\link{methyl_miss}} for creating objects,
#' \code{\link{methyl_surro}} for the input object type,
#' \code{\link{reference_fill}} and \code{\link{impute_obs}} for addressing missing data
#'
#' @name methyl_miss-class
NULL

#' imputation_stats Objects
#'
#' @description
#' The \code{imputation_stats} class stores detailed information about a missing data
#' imputation operation performed by \code{\link{impute_obs}}. These objects provide
#' comprehensive statistics about what was imputed, what was skipped, and why.
#'
#' @section Object Structure:
#' An \code{imputation_stats} object is a list with class \code{"imputation_stats"} containing:
#' \describe{
#'   \item{\code{method}}{Character string indicating the imputation method used
#'   (\code{"mean"} or \code{"median"}).}
#'   \item{\code{min_nonmiss_prop}}{Numeric value (0 to 1) indicating the minimum proportion
#'   of non-missing data required for a probe to be imputed.}
#'   \item{\code{timestamp}}{POSIXct timestamp of when imputation was performed.}
#'   \item{\code{n_total_probes}}{Total number of probes in the methylation matrix.}
#'   \item{\code{n_complete_probes}}{Number of probes with no missing values.}
#'   \item{\code{n_completely_missing_probes}}{Number of probes missing in all samples.}
#'   \item{\code{n_partially_missing_probes}}{Number of probes with some (but not all)
#'   missing values.}
#'   \item{\code{n_probes_imputed}}{Number of probes that had values imputed.}
#'   \item{\code{n_probes_skipped}}{Number of partially missing probes that were skipped
#'   because they didn't meet the \code{min_nonmiss_prop} threshold.}
#'   \item{\code{n_values_imputed}}{Total number of individual missing values that were
#'   imputed.}
#'   \item{\code{n_missing_before_imputation}}{Total missing values before imputation.}
#'   \item{\code{n_missing_after_imputation}}{Total missing values remaining after imputation.}
#'   \item{\code{imputation_rate}}{Proportion of originally missing values that were imputed
#'   (0 to 1).}
#'   \item{\code{probes_imputed}}{Character vector of CpG probe IDs that were imputed.}
#'   \item{\code{probes_skipped}}{Character vector of CpG probe IDs that were skipped.}
#'   \item{\code{values_imputed_per_probe}}{Named numeric vector showing the number of values
#'   imputed for each probe (only present when \code{return_stats = TRUE}).}
#'   \item{\code{skipped_reasons}}{Named character vector explaining why each skipped probe
#'   was not imputed (only present when \code{return_stats = TRUE}).}
#' }
#'
#' @section Creating imputation_stats Objects:
#' Objects are created automatically when calling \code{\link{impute_obs}} with
#' \code{return_stats = TRUE}:
#' \preformatted{
#' # Perform imputation with statistics
#' result <- impute_obs(
#'   my_surro,
#'   method = "mean",
#'   min_nonmiss_prop = 0.5,
#'   return_stats = TRUE
#' )
#'
#' # Access the statistics
#' stats <- result$imputation_stats
#' }
#'
#' @section Accessing Information:
#' Components can be accessed using standard list notation:
#' \preformatted{
#' # View imputation summary
#' stats$n_values_imputed
#' stats$imputation_rate
#'
#' # View which probes were imputed
#' stats$probes_imputed
#'
#' # View detailed reasons for skipped probes
#' stats$skipped_reasons
#' }
#'
#' @section Print Method:
#' The package provides a \code{print} method that displays a formatted summary:
#' \preformatted{
#' print(stats)
#' }
#' This shows an organized view of the imputation operation and its results.
#'
#' @seealso
#' \code{\link{impute_obs}} for creating objects and performing imputation,
#' \code{\link{methyl_surro}} for the input object type,
#' \code{\link{reference_fill_stats}} for reference filling statistics
#'
#' @name imputation_stats
NULL

#' reference_fill_stats Objects
#'
#' @description
#' The \code{reference_fill_stats} class stores detailed information about a reference
#' filling operation performed by \code{\link{reference_fill}}. These objects provide
#' comprehensive statistics about what was filled and how the missing data was addressed.
#'
#' @section Object Structure:
#' A \code{reference_fill_stats} object is a list with class \code{"reference_fill_stats"}
#' containing:
#' \describe{
#'   \item{\code{type}}{Character string indicating the fill type used: \code{"probes"}
#'   (only completely missing probes), \code{"obs"} (only missing observations), or
#'   \code{"all"} (both).}
#'   \item{\code{timestamp}}{POSIXct timestamp of when filling was performed.}
#'   \item{\code{n_total_probes}}{Total number of probes in the methylation matrix.}
#'   \item{\code{n_total_samples}}{Total number of samples in the methylation matrix.}
#'   \item{\code{n_completely_missing_probes_before}}{Number of completely missing probes
#'   before filling.}
#'   \item{\code{n_completely_missing_probes_after}}{Number of completely missing probes
#'   after filling.}
#'   \item{\code{n_probes_filled}}{Number of completely missing probes that were filled.}
#'   \item{\code{n_obs_filled}}{Number of individual missing observations that were filled
#'   (in partially missing probes).}
#'   \item{\code{n_total_values_filled}}{Total number of values filled (sum of probes filled
#'   × samples + individual observations filled).}
#'   \item{\code{n_missing_before}}{Total missing values before filling.}
#'   \item{\code{n_missing_after}}{Total missing values after filling.}
#'   \item{\code{fill_rate}}{Proportion of originally missing values that were filled (0 to 1).}
#'   \item{\code{probes_filled}}{Character vector of CpG probe IDs for completely missing
#'   probes that were filled.}
#'   \item{\code{obs_filled}}{Character vector of CpG probe IDs for probes that had
#'   individual observations filled.}
#'   \item{\code{values_filled_per_probe}}{Named numeric vector showing the number of values
#'   filled for each probe (only present when \code{return_stats = TRUE}).}
#'   \item{\code{probes_not_in_reference}}{Character vector of CpG probe IDs that were
#'   missing but not present in the reference values (only present when applicable).}
#' }
#'
#' @section Creating reference_fill_stats Objects:
#' Objects are created automatically when calling \code{\link{reference_fill}} with
#' \code{return_stats = TRUE}:
#' \preformatted{
#' # Perform reference filling with statistics
#' result <- reference_fill(
#'   my_surro,
#'   reference = ref_values,
#'   type = "probes",
#'   return_stats = TRUE
#' )
#'
#' # Access the statistics
#' stats <- result$reference_fill_stats
#' }
#'
#' @section Accessing Information:
#' Components can be accessed using standard list notation:
#' \preformatted{
#' # View filling summary
#' stats$n_total_values_filled
#' stats$fill_rate
#'
#' # View which probes were filled
#' stats$probes_filled
#'
#' # View detailed filling counts
#' stats$values_filled_per_probe
#' }
#'
#' @section Print Method:
#' The package provides a \code{print} method that displays a formatted summary:
#' \preformatted{
#' print(stats)
#' }
#' This shows an organized view of the reference filling operation and its results.
#'
#' @seealso
#' \code{\link{reference_fill}} for creating objects and performing reference filling,
#' \code{\link{methyl_surro}} for the input object type,
#' \code{\link{imputation_stats}} for imputation statistics
#'
#' @name reference_fill_stats
NULL
