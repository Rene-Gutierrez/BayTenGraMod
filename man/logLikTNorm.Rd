% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logLikTNorm.R
\name{logLikTNorm}
\alias{logLikTNorm}
\title{Computes the log-likelihood of the Tesor Normal Distribution}
\usage{
logLikTNorm(tensors, precisions)
}
\arguments{
\item{tensors}{A list with the sample of tensors}

\item{precisions}{A list of precision matrices.}
}
\value{
The log-likelihood for the Tensor Normal distribution precisions and
  sample specified.
}
\description{
Computes the log-likelihood of the Tesor Normal Distribution given a sample
of tensors.
}
\author{
Rene Gutierrez Marquez
}
