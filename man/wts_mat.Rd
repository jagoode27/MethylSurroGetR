% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_data.R
\name{wts_mat}
\alias{wts_mat}
\title{Weights Matrix}
\format{
A matrix with 10 rows and 3 columns.
\describe{
\item{wt_lin}{Weights for linear regression}
\item{wt_prb}{Weights for logistic (probability) regression}
\item{wt_cnt}{Weights for Poisson (count) regression}
}
}
\description{
A 10 x 3 matrix of regression weights.
}
\examples{
data(wts_mat)
str(wts_mat)
}
