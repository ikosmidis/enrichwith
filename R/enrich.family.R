#' Enrich \code{family} objects
#'
#' Enrich \code{family} objects with family- and link-specfic mathematical
#' functions
#'
#' @param object an object of class \code{\link[stats]{family}}
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with. Run \code{family_options()} to see
#'     the available options
#' @param ... currently not used
#'
#' @details
#' \code{\link{family}} objects specify the details of the models used
#' by functions such as \code{\link[stats]{glm}}. The families
#' implemented in the \code{stats} package include
#' \code{\link[stats]{binomial}}, \code{\link[stats]{gaussian}},
#' \code{\link[stats]{Gamma}}, \code{\link[stats]{inverse.gaussian}},
#' and \code{\link[stats]{poisson}}. \code{\link[stats]{family}}
#' objects specify particular characteristics of distributions from
#' the exponential family. Such distributions have probability mass or
#' density function of the form \deqn{f(y; \theta, \phi) =
#' \exp\left\{\frac{y\theta - b(\theta) + c_1(y)}{\phi/m} - a(-m/\phi)
#' + c_2(y)\right\} \quad y \in Y \subset \Re\,, \theta \in \Theta
#' \subset \Re\, , \phi > 0}{f(y, theta, phi) = exp((y * theta -
#' b(theta) + c_1(y))/(phi/m) - a(-m/phi) + c_2(y))} where \eqn{m >
#' 0}{m > 0} is an observation weight, and \eqn{a(.)}{a(.)},
#' \eqn{b(.)}{b(.)}, \eqn{c_1(.)}{c_1(.)} and \eqn{c_2(.)}{c_2(.)} are
#' sufficiently smooth, real-valued functions.
#'
#' The expected value and the variance of such distributions is
#' \eqn{\mu = b'(\theta)}{mu = b'(theta)} and \eqn{\phi V(\mu)/m}{phi * V(mu)/m},
#' repsectively, where \eqn{V(\mu)}{V(mu)} is called the variance
#' funciton. The parameter \eqn{\phi}{phi} is called a dispersion
#' parameter.
#'
#' Characteristics of the exponential family that are implementented
#' in \code{\link[stats]{family}} objects include:
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
#' The \code{enrich} can enriche \code{\link{family}} family objects
#' with extra characteristics of the family and of the chosen link
#' function. See \code{\link{enrich.link-glm}} for the enrichment of
#' link functions.
#'
#' @return The object \code{object} of class \code{\link{family}} with
#'     extra components. \code{family_options()} reutns the components
#'     and their descriptions
#' @method enrich family
#'
#' @seealso \code{\link{enrich.link-glm}}, \code{\link[stats]{family}}
#'
#' @export
#'
#' @examples
#' ## An example from ?glm to illustrate that things still work with
#' ## enriched families
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glm(counts ~ outcome + treatment, family = enrich(poisson()))
#' anova(glm.D93)
#' summary(glm.D93)
enrich.family <- function(object, with = "all", ...) {
    if (is.null(with)) {
        return(object)
    }
    what <- family_options(with)
    for (j in what) {
        ## Take care of link function specific options
        if (j %in% c("d2mu.deta", "d3mu.deta")) {
            object[[j]] <- eval(call(j, object = make.link(object$link)))
        }
        else {
            object[[j]] <- eval(call(j, object = object))
        }
    }
    object
}

#' Available options for the enrichment of \code{family} objects
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
#' family_options(what = "all")
#' family_options(print = TRUE)
#' }
#' @export
family_options <- function(what, print = missing(what)) {
    ## List the enrichment options that you would like to make
    ## avaiable for objects of class
    available_options <- c("d1variance", "d2variance", "d1afun", "d2afun", "d3afun", "d2mu.deta", "d3mu.deta",
                           "variance derivatives",
                           "function a derivatives",
                           "inverse link derivatives")
    ## Provide the descriptions of the enrichment options
    descriptions <- c("1st derivative of the variance function",
                      "2nd derivative of the variance function",
                      "1st derivative of the a function",
                      "2nd njderivative of the a function",
                      "3rd derivative of the a function",
                      "2nd derivative of the inverse link function",
                      "3rd derivative of the inverse link function",
                      "1st and 2nd derivative of the variance function",
                      "1st, 2nd and 3rd derivative of the a function",
                      "2nd and 3rd derivative of the inverse link function")
    available_options <- c(available_options, 'all')
    descriptions <- c(descriptions, 'all available options')
    if (print) {
        cat(paste(paste(available_options, descriptions, sep = " : "), "\n"))


        return(invisible())
    }
    if (any(!(what %in% available_options))) {
        stop(gettextf('one of the options %s is not implemented', paste(what, collapse =  )))
    }
    ## List what each option in available_options corresponds to
    ## (one vector of function names per option). The corresponding
    ## functions should take as input an object of class 'class'
    family_with <- list("d1variance", "d2variance", "d1afun", "d2afun", "d3afun", "d2mu.deta", "d3mu.deta",
                        c("d1variance", "d2variance"),
                        c("d1afun", "d2afun", "d3afun"),
                        c("d2mu.deta", "d3mu.deta"))
    family_with[[length(family_with) + 1]] <- unique(unlist(family_with))
    names(family_with) <- available_options
    unique(unlist(family_with[what]))
}

d1variance <- function(object) {
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

d2variance <- function(object) {
    family <- object$family
    d2variance <- switch(family,
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

d1afun <- function(object) {
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

d2afun <- function(object) {
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

d3afun <- function(object) {
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
