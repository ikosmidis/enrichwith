#' Available options for the enrichment objects of class {{class}}
#'
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
#' @examples
#' \dontrun{
#' enrichment_options.{{class}}(option = "all")
#' enrichment_options.{{class}}(all_options = TRUE)
#' }
#' @export
`get_enrichment_options.{{class}}` <- function(object, option, all_options = missing(option)) {
    ## List the enrichment options that you would like to make
    ## avaiable for objects of class
    ## enrichment_options <- {{option}}
    out <- list()
    out$option <- {{option}}
    ## Provide the descriptions of the enrichment options
    ## descriptions <- {{description}}
    out$description <- {{description}}
    ## Add all as an option
    ## enrichment_options <- c(enrichment_options, 'all')
    out$option <- c(out$option, 'all')
    ## descriptions <- c(descriptions, 'all available options')
    out$description <- c(out$description, 'all available options')
    ## components <- {{component}}
    out$component <- {{component}}
    ## components[[length(components) + 1]] <- unique(unlist(components))
    out$component[[length(out$component) + 1]] <- unique(unlist(out$component))
    ## names(components) <- enrichment_options
    names(out$component) <- out$option
    out$compute_function <- lapply(out$component, function(z) paste0('compute_', z))
    class(out) <- 'enrichment_options'
    if (all_options) {
        return(out)
    }
    ## invalid_options <- !(option %in% enrichment_options)
    invalid_options <- !(option %in% out$option)
    if (any(invalid_options)) {
        stop(gettextf('some options have not been implemented: %s', paste0('"', paste(option[invalid_options], collapse = ', '), '"')))
    }
    return(unique(unlist(out$compute_function[option])))
}


