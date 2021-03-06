% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hoff.R
\name{hoff}
\alias{hoff}
\title{hoff}
\usage{
hoff(y = y, NSCAN = 10000)
}
\arguments{
\item{y}{A numeric real-valued matrix with n number of samples and
p number of features (n > p).}

\item{NSCAN}{A numberic value decide MCMC sample size (size = NSCAN / 10)}
}
\description{
All credit to Luo, W. and Li, B. (2016),
Combining Eigenviralues and Variation of Eigenvirectors for Order
Determination, Biometrika, 103, 875–887. <doi:10.1093/biomet/asw051>
}
