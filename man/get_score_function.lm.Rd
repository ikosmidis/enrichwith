% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich.lm.R
\name{get_score_function.lm}
\alias{get_score_function.lm}
\title{Function to compute/extract a function that returns the scores
(derivatives of the log-likelihood) for an object of class
\code{lm}/\code{enriched_lm}}
\usage{
\method{get_score_function}{lm}(object, ...)
}
\arguments{
\item{object}{an object of class \code{lm} or\code{enriched_lm}}

\item{...}{currently not used}
}
\description{
Function to compute/extract a function that returns the scores
(derivatives of the log-likelihood) for an object of class
\code{lm}/\code{enriched_lm}
}
\details{
The computed/extracted function has arguments
\describe{

\item{coefficients}{the regression coefficients at which the scores
are computed. If missing then the maximum likelihood estimates are
used}

\item{dispersion}{the dispersion parameter at which the score
function is evaluated. If missing then the maximum likelihood
estimate is used}

}
}
