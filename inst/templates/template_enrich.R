#' Enrich objects of class {{class}}
#'
#' Enrich objects of class {{class}} with
#' *BRIEFLY EXPLAIN WITH WHAT*
#'
#'
#' @param object an object of class {{class}}
#' @param with a character vector of options for the enrichment of \code{object}
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
    dots <- list(...)
    enrichment_options <- get_enrichment_options(object, option = with)
    component <- unlist(enrichment_options$component)
    compute <- unlist(enrichment_options$compute_function)
    for (j in seq.int(length(component))) {
        ccall <- call(compute[j], object = object)
        for (nam in names(dots)) {
            ccall[[nam]] <- dots[[nam]]
        }
        object[[component[j]]] <- eval(ccall)
    }
    if (is.null(attr(object, "enriched"))) {
        attr(object, "enriched") <- TRUE
        classes <- class(object)
        class(object) <- c(paste0("enriched_", classes[1]), classes)
    }
    object
}



