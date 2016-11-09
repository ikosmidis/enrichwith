#' Methods to enrich list-like R objects with extra components
#'
#' The enrichwith package provides the \code{\link{enrich}} method to
#' enrich list-like R objects with new, relevant components. The
#' resulting objects preserve their class, so all methods associated
#' with them still apply. The package can also be used to produce
#' customisable source code templates for the structured
#' implementation of methods to compute new components
#'
#' Depending on the object, enriching it can be a tedious task. The
#' \pkg{enrichwith} package aims to streamline the task into 3 simple
#' steps:
#'
#' \enumerate{
#'
#' \item Use \code{\link{create_enrichwith_skeleton}} to produce a
#' customisable enrichwith template.
#'
#' \item Edit the \code{compute_*} functions by adding the specific
#' code that calculates the components.
#'
#' \item Finalise the documentation and/or include more examples.
#' }
#'
#' The first step results in a template that includes all necessary
#' functions to carry out the enrichment. The second step is where the
#' user edits the template and implements the calculation of the
#' components that the object will be enriched with. Specifically,
#' each \code{compute_*} function takes as input the object to be
#' enriched and returns the corresponding new component to be added to
#' the object.
#'
#' Everything else (for example, mapping between the enrichment
#' options and the components that the enriched object will have,
#' checks that an enrichment option exists, listing enrichment
#' options, enriching the object, and so on) is taken care of by the
#' methods in \pkg{enrichwith}.
#'
#' Developers can either put their enrichwith templates in their
#' packages or are welcome to contribute their template to enrichwith,
#' particularly if that extends core R objects.
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
#' @param with a character vector with enrichment options for \code{object}
#' @param ... Arguments to be passed to other methods
#'
#' @export
enrich <- function(object, with, ...) {
    UseMethod("enrich")
}

#' Generic method for getting available options for the enrichment
#' of objects
#'
#' @aliases print.enrichment_options
#'
#' @param object the object to be enriched
#' @param option a character vector listing the options for enriching
#'     the object
#' @param all_options if \code{TRUE} then output a data frame with the
#'     available enrichment options, their descriptions, the names of
#'     the components that each option results in, and the names of
#'     the corresponding \code{compute} functions.
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


#' Generic method for extracting or computing auxiliary functions for
#' objects
#' @param object the object to be enriched or the enriched object
#' @param ... curretly not used
get_auxiliary_functions <- function(object, ...) {
    UseMethod("get_auxiliary_functions")
}
