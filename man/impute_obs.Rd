% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_obs.R
\name{impute_obs}
\alias{impute_obs}
\title{Impute Missing Observations in methyl_surro Object}
\usage{
impute_obs(methyl_surro, method = c("mean", "median"), min_nonmiss_prop = 0)
}
\arguments{
\item{methyl_surro}{An object of class \code{methyl_surro}, containing a methylation matrix.}

\item{method}{A character string indicating the imputation method. Must be either \code{"mean"} or \code{"median"}.}

\item{min_nonmiss_prop}{The minimum proportion of non-missing data required in a probe (row) for the imputation to proceed. Must be a numeric value between 0 and 1.}
}
\value{
A \code{methyl_surro} object with missing observations imputed, updating the methylation matrix as per the given method and threshold.
}
\description{
This function imputes missing values in the methylation matrix of a \code{methyl_surro} object using either mean or median imputation.
It applies row-wise imputation based on a specified non-missing data proportion threshold.
}
\examples{
# Load the sample data
data(methyl_surro_miss)

# Apply mean imputation with a specified threshold for non-missing data in probes
result <- impute_obs(methyl_surro_miss, "mean", min_nonmiss_prop = 0.5)

# Check the imputed result
print(result$methyl)

}
