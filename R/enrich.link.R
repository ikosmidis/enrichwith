#' Enrich \code{link-glm} objects
#'
#' Enrich \code{\link[=make.link]{link-glm}} with further derivatives of
#' \code{linkinv} with respect to \code{eta}.
#'
#' @details
#'
#' @param object an object of class \code{\link{family}}
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with. Run \code{linkglm_options(print =
#'     TRUE)} to see the available options
#' @param ... currently not used
#'
#' @details
#'
#' The \code{enrich.link-glm} method supports \code{logit},
#' \code{probit}, \code{cauchit}, \code{cloglog}, \code{identity},
#' \code{log}, \code{sqrt}, \code{1/mu^2}, \code{inverse}, as well as
#' the \code{\link{power}} family of links.
#'
#' @return The object \code{object} of class
#'     \code{\link[=make.link]{link-glm}} with extra components. Run
#'     \code{linkglm_options(print = TRUE)} to see what those
#'     components are.
#'
#' @name enrich.link-glm
#' @method enrich link-glm
#' @export
#' @examples
#' elogit <- enrich(make.link("logit"), with = "inverse link derivatives")
#' str(elogit)
#' elogit$d2mu.deta
#' elogit$d3mu.deta
#'
assign(x = "enrich.link-glm",
       value = function(object, with = "all", ...) {
           if (is.null(with)) {
               return(object)
           }
           what <- linkglm_options(what = with)
           for (j in what) {
               object[[j]] <- eval(call(j, object = object))
           }
           object
       })

#' Available options for the enrichment of \code{\link[=make.link]{link-glm}} objects
#'
#' @param what a character vector listing the components that should
#'     be present in the enriched object
#' @param print if TRUE then the available enrichment options
#'     are listed
#' @details A check is being made whether the requested component is
#'     available in the \code{available_components} specification (see
#'     \code{\link{options_function}}). No check is being made on
#'     whether the functions that produce the components exist.
#' @examples
#' \dontrun{
#' linkglm_options(what = "all")
#' linkglm_options(print = TRUE)
#' }
#' @export
linkglm_options <- function(what, print = FALSE) {
    ## List the enrichment options that you would like to make
    ## avaiable for objects of class
    available_options <- c("d2mu.deta", "d3mu.deta",
                           "inverse link derivatives")
    ## Provide the descriptions of the enrichment options
    descriptions <- c("2nd derivative of the inverse link function",
                      "3rd derivative of the inverse link function",
                      "2nd and 3rd derivative of the inverse link function")
    available_options <- c(available_options, 'all')
    descriptions <- c(descriptions, 'all available options')
    if (print) {
        cat(paste(paste(available_options, descriptions, sep = ' : '), '\n'))
        return(invisible())
    }
    if (any(!(what %in% available_options))) {
        stop(gettextf('one of the options %s is not implemented', paste(what, collapse = ', ')))
	 }
    ## List what each option in available_options corresponds to
    ## (one vector of function names per option. The corresponding
    ## functions should take as input an object of class 'class'
    linkglm_with <- list("d2mu.deta", "d3mu.deta",
                         c("d2mu.deta", "d3mu.deta"),
                         c("d2mu.deta", "d3mu.deta"))
    linkglm_with[[length(linkglm_with) + 1]] <- unique(unlist(linkglm_with))
    names(linkglm_with) <- available_options
    unique(unlist(linkglm_with[what]))
}

d2mu.deta <- function(object) {
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

d3mu.deta <- function(object) {
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

if (getRversion() >= "2.15.1") globalVariables(c("lambda"))
