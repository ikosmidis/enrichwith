#' Methods to enrich various R objects with extra components
#'
#' The \pkg{enrichwith} package provides the \code{enrich} method to
#' enrich various core objects with extra, relevant components. The
#' current version has methods for enriching objects of class
#' \code{family} and \code{link-glm}. The resulting objects preserve
#' their class, so the methods associated to them still apply.
#'
#' Depending on the object, enriching it can be a tedious task. The
#' \pkg{enrichwith} package streamlines the task into 3 main steps:
#'
#' \itemize{
#'
#' \item Set up an options function. In this step the enrichment
#' options, their descriptions and the components in the enriched
#' object are specficed
#'
#' \item Author the enrichment functions. Each of those needs to have
#' the same name as the respective component in the enriched object,
#' take as input the object to be enriched, and return as output the
#' respective component for the enriched object
#'
#' \item Author the enrich method.
#'
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
enrich <- function(object, with, ...) UseMethod("enrich")
get_enrichment_options <- function(object, option, all_options) UseMethod("get_enrichment_options")

