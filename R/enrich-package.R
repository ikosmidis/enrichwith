#' Methods to enrich various R objects with extra components
#'
#' The enrichwith package provides the "enrich" method to enrich
#' various core objects with extra, relevant components. The current
#' version has methods for enriching objects of class "family" and
#' "link-glm". The resulting objects preserve their class, so the
#' methods associated to them still apply.
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
