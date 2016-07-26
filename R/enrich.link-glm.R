#' Enrich objects of class \code{\link[=make.link]{link-glm}}
#'
#'
#' Enrich objects of class \code{\link[=make.link]{link-glm}} with
#' further derivatives of \code{linkinv} with respect to \code{eta}.
#'
#' @param object an object of class \code{\link[=make.link]{link-glm}}
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with.
#' @param ... extra arguments to be passed to the \code{compute_*}
#'     functions
#'
#' @details
#' The \code{enrich.link-glm} method supports \code{logit},
#' \code{probit}, \code{cauchit}, \code{cloglog}, \code{identity},
#' \code{log}, \code{sqrt}, \code{1/mu^2}, \code{inverse}, as well as
#' the \code{\link{power}} family of links.
#'
#' @return The object \code{object} of class \code{\link[=make.link]{link-glm}}
#'     with extra components. \code{get_enrichment_options.link-glm()}
#'     returns the components and their descriptions.
#'
#' @examples
#' elogit <- enrich(make.link("logit"), with = "inverse link derivatives")
#' str(elogit)
#' elogit$d2mu.deta
#' elogit$d3mu.deta
#' @export
`enrich.link-glm` <- function(object, with = "all", ...) {
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



#' Available options for the enrichment objects of class link-glm
#'
#' @param object the object to be enriched
#' @param option a character vector listing the options for enriching
#'     the object
#' @param all_options if \code{TRUE} then output a data frame with the
#'     available enrichment options, their descriptions, the names of
#'     the components that each option results in, and the names of
#'     the corresponding \code{compute_*} functions.
#' @return an object of class \code{enrichment_options}
#'
#' @details A check is being made whether the requested option is
#'     available. No check is being made on whether the functions that
#'     produce the components exist.
#' @examples
#' \dontrun{
#' `get_enrichment_options.link-glm`(option = "all")
#' `get_enrichment_options.link-glm`(all_options = TRUE)
#' }
#' @export
`get_enrichment_options.link-glm` <- function(object, option, all_options = missing(option)) {
    ## List the enrichment options that you would like to make
    ## avaiable for objects of class
    out <- list()
    out$option <- c('d2mu.deta', 'd3mu.deta', 'inverse link derivatives')
    ## Provide the descriptions of the enrichment options
    out$description <- c('2nd derivative of the inverse link function', '3rd derivative of the inverse link function', '2nd and 3rd derivative of the inverse link function')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('d2mu.deta', 'd3mu.deta', c('d2mu.deta', 'd3mu.deta'))
    out$component[[length(out$component) + 1]] <- unique(unlist(out$component))
    names(out$component) <- out$option
    out$compute_function <- lapply(out$component, function(z) paste0('compute_', z))
    class(out) <- 'enrichment_options'
    if (all_options) {
        return(out)
    }
    invalid_options <- !(option %in% out$option)
    if (any(invalid_options)) {
        stop(gettextf('some options have not been implemented: %s', paste0('"', paste(option[invalid_options], collapse = ', '), '"')))
    }

    out <- list(option = option,
                description = out$description[out$option == option],
                component = out$component[option],
                compute_function = out$compute_function[option])
    class(out) <- 'enrichment_options'
    out
}


`compute_d2mu.deta.link-glm` <- function(object, ...) {
    mu.eta <- object$mu.eta
    linkinv <- object$linkinv
    link <- object$name
    if (grepl("mu\\^", link) & (link != "1/mu^2")) {
        d2mu.deta <- function(eta) {
            (1/lambda) * (1/lambda - 1) * eta^(1/lambda - 2)
        }
        ## set the environment so lambda can be found
        environment(d2mu.deta) <- environment(object$linkfun)
        return(d2mu.deta)
    }
    switch(link,
           "logit" = function(eta) {
               mu.eta(eta) * (1 - 2 * linkinv(eta))
           },
           "probit" = function(eta) {
               -eta * pmax(dnorm(eta),.Machine$double.eps)
           },
           "cauchit" = function(eta) {
               -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps)
           },
           "cloglog" = function(eta) {
               (1 - exp(eta)) * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
           },
           "identity" = function(eta) {
               rep.int(0, length(eta))
           },
           "log" = function(eta) {
               pmax(exp(eta), .Machine$double.eps)
           },
           "sqrt" = function(eta) {
               rep.int(2, length(eta))
           },
           "1/mu^2" = function(eta) {
               3/(4 * eta^2.5)
                                          },
           "inverse" = function(eta) {
               2/(eta^3)
           },
           ## else :
           stop(sQuote(link), " link not recognised")
           )# end switch(.)
}


`compute_d2mu.deta` <- function(object, ...) {
    UseMethod('compute_d2mu.deta')
}


`compute_d3mu.deta.link-glm` <- function(object, ...) {
    mu.eta <- object$mu.eta
    linkinv <- object$linkinv
    link <- object$name
    if (grepl("mu\\^", link) & (link != "1/mu^2")) {
        d3mu.deta <- function(eta) {
            (1/lambda) * (1/lambda - 1) * (1/lambda - 2) * eta^(1/lambda - 3)
        }
        ## set the environment so lambda can be found
        environment(d3mu.deta) <- environment(object$linkfun)
        return(d3mu.deta)
    }
    switch(link,
           "logit" = function(eta) {
               mu.eta(eta) * (1 - 6 * mu.eta(eta))
           },
           "probit" = function(eta) {
               (eta^2 - 1) * pmax(dnorm(eta),.Machine$double.eps)
           },
           "cauchit" = function(eta) {
               -2 * pi * pmax(dcauchy(eta)^2, .Machine$double.eps) + 8 * eta^2 * pi^2 * pmax(dcauchy(eta)^3, .Machine$double.eps)
           },
           "cloglog" = function(eta) {
               ((1 - exp(eta))^2 - exp(eta))  * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
           },
           "identity" = function(eta) {
               rep.int(0, length(eta))
           },
           "log" = function(eta) {
               pmax(exp(eta), .Machine$double.eps)
           },
           "sqrt" = function(eta) {
               rep.int(0, length(eta))
           },
           "1/mu^2" = function(eta) {
               -15/(8 * eta^3.5)
           },
           "inverse" = function(eta) {
               -6/(eta^4)
           },
           ## else :
           stop(sQuote(link), " link not recognised")
           )# end switch(.)
}


`compute_d3mu.deta` <- function(object, ...) {
    UseMethod('compute_d3mu.deta')
}

if (getRversion() >= "2.15.1") globalVariables(c("lambda"))






## ## Call that produced the enrichwith template for the current script:
## create_enrichwith_skeleton(class = "link-glm", option = c("d2mu.deta",
##     "d3mu.deta", "inverse link derivatives"), description = c("2nd derivative of the inverse link function",
##     "3rd derivative of the inverse link function", "2nd and 3rd derivative of the inverse link function"),
##     component = list("d2mu.deta", "d3mu.deta", c("d2mu.deta",
##         "d3mu.deta")), path = "~/Downloads", attempt_rename = TRUE)


