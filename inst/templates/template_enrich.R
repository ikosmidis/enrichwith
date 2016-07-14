#' Enrich objects of class {{class}}
#'
#' Enrich of class {{class}} with
#' *BRIEFLY EXPLAIN WITH WHAT*
#'
#'
#' @param object an object of class {{class}}
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with
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
#' @export
#' @examples
#' ## *ADD AN EXAMPLE*
`enrich.{{class}}` <- function(object, with = "all", ...) {
    if (is.null(with)) {
        return(object)
    }
    enrichment_options <- get_enrichment_options(object, option = with, ...)
    component <- unlist(enrichment_options$component)
    compute <- unlist(enrichment_options$compute_function)
    for (j in seq.int(length(component))) {
        object[[component[j]]] <- eval(call(compute[j], object = object))
    }
    object
}



