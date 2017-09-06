#' Enrich objects of class \code{\link{family}}
#'
#' Enrich objects of class \code{\link{family}} with family-specific
#' mathematical functions
#'
#' @param object an object of class \code{\link{family}}
#' @param with a character vector with enrichment options for \code{object}
#' @param ... extra arguments to be passed to the \code{compute_*}
#'     functions
#'
#' @details
#'
#' \code{\link{family}} objects specify the details of the models used
#' by functions such as \code{\link[stats]{glm}}. The families
#' implemented in the \code{stats} package include
#' \code{\link[stats]{binomial}}, \code{\link[stats]{gaussian}},
#' \code{\link[stats]{Gamma}}, \code{\link[stats]{inverse.gaussian}},
#' and \code{\link[stats]{poisson}}. \code{\link[stats]{family}}
#' objects specify particular characteristics of distributions from
#' the exponential family. Such distributions have probability mass or
#' density function of the form \deqn{f(y; \theta, \phi) =
#' \exp\left\{\frac{y\theta - b(\theta) - c_1(y)}{\phi/m} -
#' \frac{1}{2}a\left(-\frac{m}{\phi}\right) + c_2(y)\right\} \quad y
#' \in Y \subset \Re\,, \theta \in \Theta \subset \Re\, , \phi >
#' 0}{f(y, theta, phi) = exp((y * theta - b(theta) - c_1(y))/(phi/m) -
#' a(-m/phi)/2 - c_2(y))} where \eqn{m > 0}{m > 0} is an observation
#' weight, and \eqn{a(.)}{a(.)}, \eqn{b(.)}{b(.)},
#' \eqn{c_1(.)}{c_1(.)} and \eqn{c_2(.)}{c_2(.)} are sufficiently
#' smooth, real-valued functions.
#'
#' The expected value and the variance of such distributions is
#' \eqn{\mu = b'(\theta)}{mu = b'(theta)} and \eqn{\phi V(\mu)/m}{phi * V(mu)/m},
#' respectively, where \eqn{V(\mu)}{V(mu)} is called the variance
#' function. The parameter \eqn{\phi}{phi} is called a dispersion
#' parameter.
#'
#' Characteristics of the exponential family that are already
#' implemented in \code{\link[stats]{family}} objects include:
#'
#' \itemize{
#'
#' \item \code{variance}: \eqn{V(.)}{V(.)}
#'
#' \item \code{dev.resids}: \eqn{-2\left\{y c_1'(\mu) - y c_1'(y) -
#' b(c_1'(\mu)) + b(c'_1(y))\right\}}{-2(y c_1'(mu) - y c_1'(y) -
#' b(c'(mu)) + b(c'(y)))}
#'
#' \item \code{aic}: \eqn{-2\sum_{i = 1}^n \left\{\log f(y_i; \theta,
#' \phi)\right\} + 2\delta}{-2*sum(log f(y, theta, phi)) + 2*delta}
#' where \eqn{\delta}{delta} is \code{1} if the family has a dispersion
#' parameter and \code{0} else
#' }
#'
#' The \code{\link{quasi}} families differ from the other families in
#' that the variance function is not determined by the family but may
#' be supplied by the user. Also, the \code{\link{quasibinomial}} and
#' \code{\link{quasipoisson}} families differ from the
#' \code{\link{binomial}} and \code{\link{poisson}} families only in
#' that the dipsersion parameter is estimated to account for
#' overdispersion.
#'
#' For \code{\link[stats]{quasi}}, \code{dev.resid} is \eqn{m(y -
#' \mu)^2}{m*(y - mu)^2} . For \code{\link[stats]{quasibinomial}} and
#' \code{\link[stats]{quasipoisson}}, \code{dev.resid} is the same as
#' for \code{\link[stats]{quasibinomial}} and
#' \code{\link{quasipoisson}}, respectively.  The \code{aic} is
#' \code{NA} for all \code{quasi} families. See
#' \code{\link[stats]{quasi}} for more details.
#'
#' The \code{enrich} method can enrich \code{\link{family}} family
#' objects with extra characteristics of the family and of the chosen
#' link function. See \code{\link{enrich.link-glm}} for the enrichment
#' of \code{\link[=make.link]{link-glm}} objects.
#'
#' @return The object \code{object} of class \code{\link{family}} with
#'     extra components. \code{get_enrichment_options.family()}
#'     returns the components and their descriptions.
#'
#' @export
#' @examples
#'
#' ## An example from ?glm to illustrate that things still work with
#' ## enriched families
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glm(counts ~ outcome + treatment, family = enrich(poisson()))
#' anova(glm.D93)
#' summary(glm.D93)
`enrich.family` <- function(object, with = "all", ...) {
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



#' Available options for the enrichment objects of class family
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
#' get_enrichment_options.family(option = "all")
#' get_enrichment_options.family(all_options = TRUE)
#' }
#' @export
`get_enrichment_options.family` <- function(object, option, all_options = missing(option)) {
    ## List the enrichment options that you would like to make
    ## available for objects of class
    out <- list()
    out$option <- c('d1variance', 'd2variance', 'd1afun', 'd2afun', 'd3afun', 'variance derivatives', 'function a derivatives')
    ## Provide the descriptions of the enrichment options
    out$description <- c('1st derivative of the variance function', '2nd derivative of the variance function', '1st derivative of the a function', '2nd derivative of the a function', '3rd derivative of the a function', '1st and 2nd derivative of the variance function', '1st, 2nd and 3rd derivative of the a function')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('d1variance', 'd2variance', 'd1afun', 'd2afun', 'd3afun', c('d1variance', 'd2variance'), c('d1afun', 'd2afun', 'd3afun'))
    out$component[[length(out$component) + 1]] <- unique(unlist(out$component))
    names(out$component) <- names(out$description) <- out$option
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
                description = out$description[option],
                component = out$component[option],
                compute_function = out$compute_function[option])
    class(out) <- 'enrichment_options'
    out
}


`compute_d1variance.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "poisson" = function(mu) {
               rep(1, length(mu))
           },
           "quasipoisson" = function(mu) {
               rep(1, length(mu))
           },
           "gaussian" = function(mu) {
               rep(0, length(mu))
                                },
           "binomial" = function(mu) {
               1 - 2 * mu
                     },
           "quasibinomial" = function(mu) {
               1 - 2 * mu
           },
           "Gamma" = function(mu) {
               2 * mu
           },
           "inverse.gaussian" = function(mu) {
               3 * mu^2
           },
           "quasi" = switch(object$varfun,
                            "constant" = function(mu) {
                                rep(0, length(mu))
                            },
                            "mu(1-mu)" = function(mu) {
                                1 - 2 * mu
                            },
                            "mu" = function(mu) {
                                rep(1, length(mu))
                            },
                            "mu^2" = function(mu) {
                                2 * mu
                            },
                            "mu^3" = function(mu) {
                                3 * mu^2
                            },
                            stop(sQuote(object$varfun), " variance function not supported")),
           stop(sQuote(family), " family not supported"))
}


`compute_d1variance` <- function(object, ...) {
    UseMethod('compute_d1variance')
}


`compute_d2variance.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "poisson" = function(mu) {
               rep(0, length(mu))
           },
           "quasipoisson" = function(mu) {
               rep(0, length(mu))
           },
           "gaussian" = function(mu) {
                             rep(0, length(mu))
           },
           "binomial" = function(mu) {
               rep(-2, length(mu))
           },
           "quasibinomial" = function(mu) {
               rep(-2, length(mu))
           },
           "Gamma" = function(mu) {
               rep(2, length(mu))
           },
           "inverse.gaussian" = function(mu) {
               6 * mu
           },
           "quasi" = switch(object$varfun,
                            "constant" = function(mu) {
                                rep(0, length(mu))
                            },
                            "mu(1-mu)" = function(mu) {
                                rep(-2, length(mu))
                            },
                            "mu" = function(mu) {
                                rep(0, length(mu))
                                          },
                            "mu^2" = function(mu) {
                                rep(2, length(mu))
                            },
                            "mu^3" = function(mu) {
                                6 * mu
                            },
                            stop(sQuote(object$varfun), " variance function not supported")),
           stop(sQuote(family), " family not supported"))
}


`compute_d2variance` <- function(object, ...) {
    UseMethod('compute_d2variance')
}


`compute_d1afun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "gaussian" = function(zeta) {
               -1/zeta
           },
           "Gamma" = function(zeta) {
               -2*psigamma(-zeta, 0) + 2*log(-zeta)
           },
           "inverse.gaussian" = function(zeta) {
               -1/zeta
           })
}


`compute_d1afun` <- function(object, ...) {
    UseMethod('compute_d1afun')
}


`compute_d2afun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "gaussian" = function(zeta) {
               1/zeta^2
           },
           "Gamma" = function(zeta) {
               2*psigamma(-zeta, 1) + 2/zeta
           },
           "inverse.gaussian" = function(zeta) {
               1/zeta^2
           })
}


`compute_d2afun` <- function(object, ...) {
    UseMethod('compute_d2afun')
}


`compute_d3afun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "gaussian" = function(zeta) {
               -2/zeta^3
           },
           "Gamma" = function(zeta) {
               -2*psigamma(-zeta, 2) - 2/zeta^2
           },
           "inverse.gaussian" = function(zeta) {
               -2/zeta^3
           })
}


`compute_d3afun` <- function(object, ...) {
    UseMethod('compute_d3afun')
}






## ## Call that produced the enrichwith template for the current script:
## create_enrichwith_skeleton(class = "family", option = c("d1variance",
##     "d2variance", "d1afun", "d2afun", "d3afun", "variance derivatives",
##     "function a derivatives"), description = c("1st derivative of the variance function",
##     "2nd derivative of the variance function", "1st derivative of the a function",
##     "2nd derivative of the a function", "3rd derivative of the a function",
##     "1st and 2nd derivative of the variance function", "1st, 2nd and 3rd derivative of the a function"),
##     component = list("d1variance", "d2variance", "d1afun", "d2afun",
##         "d3afun", c("d1variance", "d2variance"), c("d1afun",
##             "d2afun", "d3afun")), path = "~/Downloads", attempt_rename = TRUE)



