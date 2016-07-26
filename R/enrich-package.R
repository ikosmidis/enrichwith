#' Methods to enrich various R objects with extra components
#'
#' The \pkg{enrichwith} package provides the \code{enrich} method to
#' enrich various core objects with extra, relevant components. The
#' current version has methods for enriching objects of class
#' \code{family} and \code{link-glm}. The resulting objects preserve
#' their class, so all methods associated to them still apply.
#'
#' Depending on the object, enriching it can be a tedious task. The
#' \pkg{enrichwith} package is an attempt to structure the task into 3 simple steps:
#'
#' \enumerate{
#'
#' \item Use \code{\link{create_enrichwith_skeleton}} to produce an
#' enrichwith template.
#'
#' \item Complete the functions for the \code{compute_*} methods in the
#' enrichwith template.
#'
#' \item Finalise the documentation and/or include more examples.
#' }
#'
#' The first step results in a template that includes all necessary
#' functions to carry out the enrichment. The second step is where the
#' user writes the code to calculate the components that the object
#' will be enriched with. Specifically, each \code{compute_*} function
#' takes as input the object to be enriched and return the
#' corresponding new component to be used for the enrichment of the
#' object. Everything else (for example, mapping between the
#' enrichment options and the components that the enriched object will
#' have, checks that an enrichment option exists, listing enrichment
#' options, enriching the object, and so on) is taken care of by the
#' methods in \pkg{enrichwith}.
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

