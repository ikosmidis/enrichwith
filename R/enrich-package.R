#' Methods to enrich various R objects with extra components
#'
#' The \pkg{enrichwith} package provides the \code{enrich} method to
#' enrich various core objects with extra, relevant components. The
#' current version has methods for enriching objects of class
#' \code{family} and \code{link-glm}. The resulting objects preserve
#' their class, so the methods associated to them still apply.
#'
#' Depending on the object, enriching it can be a tedious task. The
#' \pkg{enrichwith} package streamlines the task into 3 simple steps:
#'
#' \itemize{
#'
#' \item Use \code{\link{create_enrichwith_skeleton}} to produce an
#' enrichwith template. This template includes all necessary functions
#' to carry out the enrichment.
#'
#' \item Write the appropriate code for the \code{compute_*} methods
#'
#' \item Finalise the documentation and/or include more examples
#' }
#'
#' @docType package
#' @name enrichwith
#' @import stats
#'
NULL
#> NULL

## register S3 methods
#' Generic method for enriching objects
#'
#' @param object the object to be enriched
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with
#' @param ... Arguments to be passed to other methods
#'
#' @export
enrich <- function(object, with, ...) {
    UseMethod("enrich")
}

#' Generic method for getting available options for the enrichment
#' objects
#'
#' @aliases print.enrichment_options
#'
#' @param object the object to be enriched
#' @param option a character vector listing the options for enriching
#'     the object
#' @param all_options if \code{TRUE} then output a data frame with the
#'     available enrichment options, their descriptions, the names of
#'     the components that each option results in, and the names of
#'     the corresponding \code{compute} funcitons.
#' @return if \code{all_options = TRUE} then an object of class
#'     \code{enrichment_options} is returned, otherwise if
#'     \code{option} is specified the output is a character vector
#'     with the names of the functions that compute the enrichment
#'     components
#'
#' @details A check is being made whether the requested option is
#'     available. No check is being made on whether the functions that
#'     produce the components exist.
#' @export
get_enrichment_options <- function(object, option, all_options) {
    UseMethod("get_enrichment_options")
}

