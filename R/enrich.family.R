#' Enrich \code{family} objects
#'
#' Enrich \code{family} objects with family- and link-specfic mathematical
#' functions
#'
#' @param object an object of class \code{\link{family}}
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with. Run \code{family_options(print =
#'     TRUE)} to see the available options
#' @param ... currently not used
#'
#' @details \code{\link{family}} objects specify the details of the
#'     models used by functions such as \code{\link{glm}}. The
#'     families implemented in the stats package include
#'     \code{\link{binomial}}, \code{\link{poisson}},
#'     \code{\link{Gamma}}, \code{\link{inverse.gaussian}},
#'     \code{\link{gaussian}}, \code{\link{quasi}},
#'     \code{\link{quasibinomial}}, \code{\link{quasipoisson}}. All
#'     but the quasi families attempt to specify certain
#'     characteristics of distributions from the exponential family,
#'     which have probability mass or density function of the form
#' \deqn{f(y; \theta) = \exp\left\{\frac{y\theta - b(\theta) + c_1(y)}{phi/m} - a(-m/\phi) + c_2(y)\right\}}{}
#' where \eqn{m} is an observation weight.
#'
#' These families correspond to models with \eqn{\mu = E(Y; \theta) =
#' b'(\theta)} and \eqn{Var(Y; \theta) = \phi V(\mu)/m} (see
#' \code{\link{family}} where \code{variance} corresponds to
#' \eqn{V(\mu)} above).
#'
#' Then the deviance residuals (see \code{dev.resids}) in
#' \code{\link{family}}) are
#' \deqn{-2\left\{y c_1'(\mu) - y c_1'(y) - b(c'(\mu)) + b(c'(y))\right\}}
#' with the second and third derivatives of
#' the inverse link function with respect to \code{eta}
#' (\code{d2mu.eta} and \code{d3mu.eta}; the first is \code{mu.eta}),
#' the first and second derivative of the variance function with
#' respect to \code{mu} (\code{d1variance}, \code{d2variance}), and
#' the first three derivatives of
#'
#' @return The object \code{object} of class \code{\link{family}} with
#'     extra components. Run \code{family_options(print = TRUE)} to
#'     see what those components are.
#' @method enrich family
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
family_options <- function(what, print = FALSE) {
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
                      "2nd derivative of the a function",
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
