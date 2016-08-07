#' Enrich objects of class glm
#'
#' Enrich of class glm with
#' *BRIEFLY EXPLAIN WITH WHAT*
#'
#'
#' @param object an object of class glm
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with
#' @param ... extra arguments to be passed to the
#'     \code{compute_*} functions
#'
#' @details
#' *ADD DETAILS IF NECESSARY*
#'
#' @return The object \code{object} of class glm with extra
#'     components. \code{get_enrichment_options.glm()} returns
#'     the components and their descriptions.
#'
#' @export
#' @examples
#' ## *ADD AN EXAMPLE*
`enrich.glm` <- function(object, with = "all", ...) {
    if (is.null(with)) {
        return(object)
    }
    enrichment_options <- get_enrichment_options(object, option = with, ...)
    component <- unlist(enrichment_options$component)
    compute <- unlist(enrichment_options$compute_function)
    for (j in seq.int(length(component))) {
        object[[component[j]]] <- eval(call(compute[j], object = object))
    }
    if (is.null(attr(object, "enriched"))) {
        attr(object, "enriched") <- TRUE
        classes <- class(object)
        class(object) <- c(paste0("enriched_", classes[1]), classes)
    }
    object
}


#' Available options for the enrichment objects of class glm
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
#' get_enrichment_options.glm(option = "all")
#' get_enrichment_options.glm(all_options = TRUE)
#' }
#' @export
`get_enrichment_options.glm` <- function(object, option, all_options = missing(option)) {
    ## List the enrichment options that you would like to make
    ## avaiable for objects of class
    out <- list()
    out$option <- c('first-order bias', 'observed information', 'MLE of dispersion')
    ## Provide the descriptions of the enrichment options
    out$description <- c('the first term in the expansion of the bias of the maximum likelihood estimator', 'the observed information matrix evaluated at the maximum likelihood estimates', 'the maximum likelihood estimator of the dispersion parameter')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('bias', 'observed.information', 'dispersion.mle')
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


`compute_bias.glm` <- function(object, ...) {
    ## Enrich object with the mle of diseprsion
    object <- enrich(object, with = "MLE of dispersion")
    ## Enrich link-glm object with the 2nd and 3rd derivatives of the inverse link function
    link <- enrich(make.link(object$family$link), with = "inverse link derivatives")
    if (family$family %in% c("poisson", "binomial")) {
        dispersionML <- 1
    }
    else {
        y <- model.response(object$model)
        mus <- fitted(object)
        weights <- weights(object, type = "prior")
        nobs <- length(mus)
        keep <- weights > 0
        dfResidual <- sum(keep) - object$qr$rank
    }

    X <- model.matrix(object)
    W <- mode.
}


`compute_bias` <- function(object, ...) {
    UseMethod('compute_bias')
}


`compute_observed.information.glm` <- function(object, ...) {
    ## Write some code to compute the component observed.information using
    ## the components of object and any other arguments in ...
    ## For example
    cat('component', 'observed.information', 'for objects of class', 'glm', "\n")
}


`compute_observed.information` <- function(object, ...) {
    UseMethod('compute_observed.information')
}


`compute_dispersion.mle.glm` <- function(object, ...) {
    if (object$method != "glm.fit") {
        stop("method is not 'glm.fit'")
    }
    if (is.null(object$model)) {
        object <- update(object, model = TRUE)
    }
    ## Enrich family object
    family <- enrich(object$family, with = "all")
    d1afun <- family$d1afun
    d2afun <- family$d2afun
    d3afun <- family$d3afun
    if (family$family %in% c("poisson", "binomial")) {
        dispersionML <- 1
    }
    else {
        y <- model.response(object$model)
        mus <- fitted(object)
        weights <- weights(object, type = "prior")
        nobs <- length(mus)
        keep <- weights > 0
        dfResidual <- sum(keep) - object$qr$rank
        gradfun <- function(disp) {
            prec <- 1/disp
            zetas <- -weights * prec
            ## Evaluate the derivatives of the a function only for
            ## observations with non-zero weight
            d1afuns <- d2afuns <- d3afuns <- rep(NA, nobs)
            d1afuns[keep] <- d1afun(zetas[keep])
            d2afuns[keep] <- d2afun(zetas[keep])
            d3afuns[keep] <- d3afun(zetas[keep])
            devianceResiduals <- family$dev.resids(y, mus, weights)
            EdevianceResiduals <- weights * d1afuns
            prec^2 * sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / 2
        }
        if (dfResidual > 0) {
            dispFit <- try(uniroot(f = gradfun, lower = .Machine$double.eps, upper = 10000, tol = 1e-06), silent = FALSE)
            if (inherits(dispFit, "try-error")) {
                warning("the ML estimate of the dispersion could not be calculated")
                dispersionML <- NA
            }
            else {
                dispersionML <- dispFit$root
            }
        }
        else { ## if the model is saturated dispersionML is NA
            dispersionML <- NA
        }
    }
    dispersionML
}


`compute_dispersion.mle` <- function(object, ...) {
    UseMethod('compute_dispersion.mle')
}






## ## Call that produced the enrichwith template for the current script:
## create_enrichwith_skeleton(class = "glm", option = c("first-order bias",
##     "observed information", "MLE of dispersion"), description = c("the first term in the expansion of the bias of the maximum likelihood estimator",
##     "the observed information matrix evaluated at the maximum likelihood estimates",
##     "the maximum likelihood estimator of the dispersion parameter"),
##     component = list("bias", "observed.information", "dispersion.mle"),
##     path = "~/Downloads", attempt_rename = TRUE)

