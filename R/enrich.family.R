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
#' \code{\link[stats]{family}} objects specify characteristics of the
#' models used by functions such as \code{\link[stats]{glm}}. The
#' families implemented in the \code{stats} package include
#' \code{\link[stats]{binomial}}, \code{\link[stats]{gaussian}},
#' \code{\link[stats]{Gamma}}, \code{\link[stats]{inverse.gaussian}},
#' and \code{\link[stats]{poisson}}, which are all special cases of
#' the exponential family of distributions that have probability mass
#' or density function of the form \deqn{f(y; \theta, \phi) =
#' \exp\left\{\frac{y\theta - b(\theta) - c_1(y)}{\phi/m} -
#' \frac{1}{2}a\left(-\frac{m}{\phi}\right) + c_2(y)\right\} \quad y
#' \in Y \subset \Re\,, \theta \in \Theta \subset \Re\, , \phi >
#' 0}{f(y, theta, phi) = exp((y * theta - b(theta) - c_1(y))/(phi/m) -
#' a(-m/phi)/2 + c_2(y))} where \eqn{m > 0}{m > 0} is an observation
#' weight, and \eqn{a(.)}{a(.)}, \eqn{b(.)}{b(.)},
#' \eqn{c_1(.)}{c_1(.)} and \eqn{c_2(.)}{c_2(.)} are sufficiently
#' smooth, real-valued functions.
#'
#' The current implementation of \code{\link[stats]{family}} objects
#' includes the variance function (\code{variance}), the deviance
#' residuals (\code{dev.resids}), and the Akaike information criterion
#' (\code{aic}). See, also \code{\link{family}}.
#'
#' The \code{enrich} method can further enrich exponential
#' \code{\link{family}} distributions with \eqn{\theta}{theta} in
#' terms of \eqn{\mu}{mu} (\code{theta}), the functions
#' \eqn{b(\theta)}{b(theta)} (\code{bfun}), \eqn{c_1(y)}{c_1(y)}
#' (\code{c1fun}), \eqn{c_2(y)}{c_2(y)} (\code{c2fun}),
#' \eqn{a(\zeta)}{a(zeta)} (\code{fun}), the first two derivatives of
#' \eqn{V(\mu)}{V(mu)} (\code{d1variance} and \code{d2variance},
#' respectively), and the first four derivatives of
#' \eqn{a(\zeta)}{a(zeta)} (\code{d1afun}, \code{d2afun},
#' \code{d3afun}, \code{d4afun}, respectively).
#'
#' Corresponding enrichment options are also avaialble for
#' \code{\link[stats]{quasibinomial}},
#' \code{\link[stats]{quasipoisson}} and \code{\link[gnm]{wedderburn}}
#' families.
#'
#' The \code{\link[stats]{quasi}} families are enriched with
#' \code{d1variance} and \code{d2variance}.
#'
#' See \code{\link{enrich.link-glm}} for the enrichment of
#' \code{\link[=make.link]{link-glm}} objects.
#'
#' @return The object \code{object} of class \code{\link{family}} with
#'     extra components. \code{get_enrichment_options.family()}
#'     returns the components and their descriptions.
#'
#' @seealso \code{\link{enrich.link-glm}}, \code{\link{make.link}}
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
    out$option <- c('theta', 'bfun', 'c1fun', 'c2fun', 'd1variance', 'd2variance', 'afun', 'd1afun', 'd2afun', 'd3afun', 'd4fun', 'variance derivatives', 'function a derivatives')
    ## Provide the descriptions of the enrichment options
    out$description <- c('natural parameter', 'cumulant transform', 'c1 function', 'c2 function', '1st derivative of the variance function', '2nd derivative of the variance function', 'a function', '1st derivative of the a function', '2nd derivative of the a function', '3rd derivative of the a function', '4th derivative of the a function', '1st and 2nd derivative of the variance function', '1st, 2nd and 3rd derivative of the a function')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('theta', 'bfun', 'c1fun', 'c2fun', 'd1variance', 'd2variance', 'afun', 'd1afun', 'd2afun', 'd3afun', 'd4afun', c('d1variance', 'd2variance'), c('d1afun', 'd2afun', 'd3afun'))
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


`compute_theta.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "poisson" = function(mu) {
               log(mu)
           },
           "quasipoisson" = function(mu) {
               log(mu)
           },
           "gaussian" = function(mu) {
               mu
           },
           "binomial" = function(mu) {
               log(mu/(1 - mu))
           },
           "quasibinomial" = function(mu) {
               log(mu/(1 - mu))
           },
           "wedderburn" = function(mu) {
               log(mu/(1 - mu))
           },
           "Gamma" = function(mu) {
               -1/mu
           },
           "inverse.gaussian" = function(mu) {
               -1/2*mu^2
           })
}

`compute_theta` <- function(object, ...) {
    UseMethod('compute_theta')
}


`compute_bfun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "poisson" = function(theta) {
               exp(theta)
           },
           "quasipoisson" = function(theta) {
               exp(theta)
           },
           "gaussian" = function(theta) {
               theta^2/2
           },
           "binomial" = function(theta) {
               log(1 + exp(theta))
           },
           "quasibinomial" = function(theta) {
               log(1 + exp(theta))
           },
           "wedderburn" = function(theta) {
               log(1 + exp(theta))
           },
           "Gamma" = function(theta) {
               -log(-theta)
           },
           "inverse.gaussian" = function(theta) {
               -sqrt(-2*theta)
           })
}

`compute_bfun` <- function(object, ...) {
    UseMethod('compute_bfun')
}



`compute_c1fun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "poisson" = function(y) {
        0
    },
    "quasipoisson" = function(y) {
        0
    },
    "gaussian" = function(y) {
        y^2/2
    },
    "binomial" = function(y) {
        0
    },
    "quasibinomial" = function(y) {
        0
    },
    "wedderburn" = function(y) {
        0
    },
    "Gamma" = function(y) {
        -log(y)
    },
    "inverse.gaussian" = function(y) {
        1/(2*y)
    })
}

`compute_c1fun` <- function(object, ...) {
    UseMethod('compute_c1fun')
}



`compute_c2fun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "poisson" = function(y) {
        -log(gamma(y+1))
    },
    "quasipoisson" = function(y) {
        -log(gamma(y+1))
    },
    "gaussian" = function(y) {
        -log(2*pi)/2
    },
    "binomial" = function(y, m) {
        log(choose(m, m * y))
    },
    "quasibinomial" = function(y, m) {
        log(choose(m, m * y))
    },
    "wedderburn" = function(y, m) {
        log(choose(m, m * y))
    },
    "Gamma" = function(y) {
        -log(y)
    },
    "inverse.gaussian" = function(y) {
        -0.5 * log(2 * pi * y^3)
    })
}

`compute_c2fun` <- function(object, ...) {
    UseMethod('compute_c2fun')
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
           "wedderburn" = function(mu) {
               2 * mu * (1 - mu) * (1 - 2 * mu)
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
           "wedderburn" = function(mu) {
               2 - 12 * mu * ( 1- mu)
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

`compute_afun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "gaussian" = function(zeta) {
        -log(-zeta)
    },
    "Gamma" = function(zeta) {
        2*log(gamma(-zeta)) + 2*zeta*log(-zeta)
    },
    "inverse.gaussian" = function(zeta) {
        -log(-zeta)
    },
    "poisson" = function(zeta) {
        0
    },
    "quasipoisson" = function(zeta) {
        0
    },
    "binomial" = function(zeta) {
        0
    },
    "quasibinomial" = function(zeta) {
        0
    },
    "wedderburn" = function(zeta) {
        0
    })
}


`compute_afun` <- function(object, ...) {
    UseMethod('compute_afun')
}


`compute_d1afun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "gaussian" = function(zeta) {
        -1/zeta
    },
    "Gamma" = function(zeta) {
        -2*psigamma(-zeta, 0) + 2*log(-zeta) + 2
        ## this is the expectation of dev.resids + 2, because of the way dev.resids is implemented
    },
    "inverse.gaussian" = function(zeta) {
        -1/zeta
    },
    "poisson" = function(zeta) {
        0
    },
    "quasipoisson" = function(zeta) {
        0
    },
    "binomial" = function(zeta) {
        0
    },
    "quasibinomial" = function(zeta) {
        0
    },
    "wedderburn" = function(zeta) {
        0
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
    },
    "poisson" = function(zeta) {
        0
    },
    "quasipoisson" = function(zeta) {
        0
    },
    "binomial" = function(zeta) {
        0
    },
    "quasibinomial" = function(zeta) {
        0
    },
    "wedderburn" = function(zeta) {
        0
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
    },
    "poisson" = function(zeta) {
        0
    },
    "quasipoisson" = function(zeta) {
        0
    },
    "binomial" = function(zeta) {
        0
    },
    "quasibinomial" = function(zeta) {
        0
    },
    "wedderburn" = function(zeta) {
        0
    })
}

`compute_d3afun` <- function(object, ...) {
    UseMethod('compute_d3afun')
}


`compute_d4afun.family` <- function(object, ...) {
    family <- object$family
    switch(family,
           "gaussian" = function(zeta) {
        6/zeta^4
    },
    "Gamma" = function(zeta) {
        2 * psigamma(-zeta, 3) + 4/zeta^3
    },
    "inverse.gaussian" = function(zeta) {
        6/zeta^4
    },
    "poisson" = function(zeta) {
        0
    },
    "quasipoisson" = function(zeta) {
        0
    },
    "binomial" = function(zeta) {
        0
    },
    "quasibinomial" = function(zeta) {
        0
    },
    "wedderburn" = function(zeta) {
        0
    })
}


`compute_d4afun` <- function(object, ...) {
    UseMethod('compute_d4afun')
}

## ## Call that produced the initial enrichwith template for the current script:
## create_enrichwith_skeleton(class = "family", option = c("d1variance",
##     "d2variance", "d1afun", "d2afun", "d3afun", "variance derivatives",
##     "function a derivatives"), description = c("1st derivative of the variance function",
##     "2nd derivative of the variance function", "1st derivative of the a function",
##     "2nd derivative of the a function", "3rd derivative of the a function",
##     "1st and 2nd derivative of the variance function", "1st, 2nd and 3rd derivative of the a function"),
##     component = list("d1variance", "d2variance", "d1afun", "d2afun",
##         "d3afun", c("d1variance", "d2variance"), c("d1afun",
##             "d2afun", "d3afun")), path = "~/Downloads", attempt_rename = TRUE)



