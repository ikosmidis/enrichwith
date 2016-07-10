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
#' \item Run \code{\link{build_enrichwith_skeleton}}. This step will
#' create a file where all functions for an \code{\link{enrich}}
#' template are set-up
#'
#' \item Write the \code{compute_*} methods
#'
#' \item Complete the documentation and/or examples
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

