#' Enrich objects of class {{class}}
#'
#' Enrich of class {{class}} with
#' *BRIEFLY EXPLAIN WITH WHAT*
#'
#' @details
#'
#' @param object an object of class {{class}}
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with.
#' @param ... extra arguments to be passed to the
#'     \code{compute_*} functions
#'
#' @details
#' *ADD DETAILS IF NECESSARY*
#'
#' @return The object \code{object} of class {{class}} with extra
#'     components. \code{get_enrichment_options.{{class}}()} returns
#'     the components and their descriptions.
#'
#' @name enrich.{{class}}
#' @method enrich {{class}}
#' @export
#' @examples
#' ## *ADD AN EXAMPLE*
`enrich.{{class}}` <- function(object, with = "all", ...) {
           if (is.null(with)) {
               return(object)
           }
           what <- get_enrichment_options(object, option = with, ...)
           for (j in what) {
               object[[j]] <- eval(call(j, object = object))
           }
           object
}



