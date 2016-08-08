#' Enrich objects of class \code{\link{glm}}
#'
#' Enrich objects of class \code{\link{glm}} with any or all of a set
#' auxiliary functions, the maximum likelihood estimate of the
#' dispersion parameter, the expected or observed information at the
#' maximum likelihood estimator, and the first term in the expansion
#' of the bias of the maximum likelihood estimator.
#'
#' @param object an object of class glm
#' @param with a character vector with the names of the components to
#'     enrich \code{object} with
#' @param ... extra arguments to be passed to the
#'     \code{compute_*} functions
#'
#' @details
#'
#' The auxiliary functions consist of the score functions, the
#' expected or observed information and the first-order bias of the
#' maximum likelihood estimator as functions of the model parameters.
#'
#' @return
#'
#' The object \code{object} of class \code{\link{glm}} with extra
#' components. \code{get_enrichment_options.glm()} returns the
#' components and their descriptions.
#'
#' @export
#' @examples
#' \dontrun{
#'
#' ## Reproduce left plot in Figure 4.1 in Kosimdis (2007)
#' ## (see http://www.ucl.ac.uk/~ucakiko/files/ikosmidis_thesis.pdf)
#' mod <- glm(1 ~ 1, weights = 10, family = binomial())
#' enriched_mod <- enrich(aa, with = "auxiliary functions")
#' biasfun <- enriched_mod$auxiliary_functions$bias
#' probabilities <- seq(1e-02, 1 - 1e-02, length = 100)
#' biases <- Vectorize(biasfun)(qlogis(probabilities))[1,]
#' plot(probabilities, biases, type = "l", ylim = c(-0.5, 0.5)
#'      xlab = expression(pi), ylab = "first-order bias")
#' abline(h = 0, lty = 2)
#' title("First-order bias of the MLE of the log-odds", sub = "m = 10")
#' }
`enrich.glm` <- function(object, with = "all", ...) {
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



#' Available options for the enrichment objects of class
#' \code{\link{glm}}
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
    ## available for objects of class
    out <- list()
    out$option <- c('auxiliary functions', 'score vector', 'mle of dispersion', 'expected information', 'observed information', 'first-order bias')
    ## Provide the descriptions of the enrichment options
    out$description <- c('various likelihood-based quantities as functions of the model parameters', 'gradient of the log-likelihood at the mle', 'mle of the dispersion parameter', 'expected information matrix evaluated at the mle', 'observed information matrix evaluated at the mle', 'first term in the expansion of the bias of the mle at the mle')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('auxiliary_functions', 'score_mle', 'dispersion_mle', 'expected_information_mle', 'observed_information_mle', 'bias_mle')
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


`compute_auxiliary_functions.glm` <- function(object, ...) {
    if (is.null(object$model)) {
        object <- update(object, model = TRUE)
    }

    ## Enrich link-glm and family objects
    link <- enrich(make.link(object$family$link), with = "all")
    family <- enrich(object$family, with = "all")

    ## Extract functions from link
    linkfun <- link$linkfun
    linkinv <- link$linkinv
    d1mu <- link$mu.eta
    d2mu <- link$d2mu.deta

    ## Extract functions from family
    variance <- family$variance
    d1variance <- family$d1variance
    d1afun <- family$d1afun
    d2afun <- family$d2afun
    d3afun <- family$d3afun
    dev.resids <- family$dev.resids
    aic <- family$aic

    ## Get x, y, ...
    y <- object$y ## We do not do model.response here to cover the
                  ## case where successes and failures are give for a
                  ## binomial glm
    x <- model.matrix(object)
    nobs <- nobs(object)
    nvar <- ncol(x)
    off <- model.offset(object$model)
    prior_weights <- weights(object, type = "prior")
    keep <- prior_weights > 0
    dfResidual <- sum(keep) - object$rank

    if (is.null(off)) {
        off <- rep(0, nobs)
    }

    score <- function(coefs, dispersion = 1) {
        predictors <- drop(x %*% coefs + off)
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        ## Score for coefss
        score_coefs <- colSums(prior_weights * d1mus * (y - fitted_values) * x / variances)/dispersion
        ## Score for dispersion
        if (family$family %in% c("poisson", "binomial")) {
            score_dispersion <- 0
        }
        else {
            zetas <- -prior_weights/dispersion
            d1afuns <- rep(NA, nobs)
            d1afuns[keep] <- d1afun(zetas[keep])
            devianceResiduals <- family$dev.resids(y, fitted_values, prior_weights)
            EdevianceResiduals <- prior_weights * d1afuns
            score_dispersion <- sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / (2 * dispersion^2)
        }
        ## Overall score
        out <- c(score_coefs, score_dispersion)
        names(out) <- paste0("grad_", c(names(score_coefs), "dispersion"))
        out
    }

    information <- function(coefs, dispersion = 1,
                            type = c("expected", "observed"), QR = FALSE) {
        type <- match.arg(type)
        predictors <- drop(x %*% coefs + off)
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        working_weights <- prior_weights * d1mus^2 / variances
        wx <- x * sqrt(working_weights)
        if (QR) {
            return(qr(wx))
        }
        ## expected info coefs-coefs
        info_coefs <- crossprod(wx) / dispersion
        if (type == "observed") {
            d2mus <- d2mu(predictors)
            d1variances <- d1variance(fitted_values)
            w <- prior_weights * (d2mus / variances - d1mus^2 * d1variances / variances^2) * (y - fitted_values)
            ## observed info coefs-coefs
            info_coefs <- info_coefs - t(x * w) %*% x / dispersion
        }
        rownames(info_coefs) <- colnames(info_coefs) <- colnames(x)
        ## If there is no dispersion parameter then return the
        ## information for the coefficients only
        if (family$family %in% c("poisson", "binomial")) {
            return(info_coefs)
        }
        ## If there is a dispersion parameter then return the
        ## information on the coefficients and tha dispersion
        else {
            ## expected info coefs-dispersion
            info_cross <- rep(0, ncol(info_coefs))
            ## expected info dispersion-dispersion
            zetas <- -prior_weights/dispersion
            d2afuns <- rep(NA, nobs)
            d2afuns[keep] <- d2afun(zetas[keep])
            info_dispe <- sum(prior_weights^2 * d2afuns, na.rm = TRUE) / (2 * dispersion^4)
            if (type == "observed") {
                ## observed info coefs-dispersion
                info_cross <- info_cross + colSums(prior_weights * d1mus * (y - fitted_values) * x / variances)/dispersion^2
                ## observed info dispersion-dispersion
                zetas <- -prior_weights/dispersion
                d1afuns <- rep(NA, nobs)
                d1afuns[keep] <- d1afun(zetas[keep])
                devianceResiduals <- family$dev.resids(y, fitted_values, prior_weights)
                EdevianceResiduals <- prior_weights * d1afuns
                info_dispe <- info_dispe + sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / dispersion^3
            }
            out <- rbind(cbind(info_coefs, info_cross),
                         c(info_cross, info_dispe))
            colnames(out) <- rownames(out) <- c(colnames(x), "dispersion")
            out
        }
    }

    bias <- function(coefs, dispersion = 1) {
        predictors <- drop(x %*% coefs + off)
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        d2mus <- d2mu(predictors)
        variances <- variance(fitted_values)
        working_weights <- prior_weights * d1mus^2 / variances
        Qr <- information(coefs, dispersion = dispersion, QR = TRUE)
        Q <- qr.Q(Qr)
        hats <- rowSums(Q * Q)
        ksi <- -0.5 * dispersion * d2mus * hats / (d1mus * sqrt(working_weights))
        bias_beta <- drop(tcrossprod(ksi %*% Q, solve(qr.R(Qr))))
        if (family$family %in% c("poisson", "binomial")) {
            bias_dispersion <- 0
        }
        else {
            if (dfResidual > 0) {
                ## Enrich family object with the the derivatives of the a
                ## function (see ?enrich.family for details)
                zetas <- -prior_weights/dispersion
                d2afuns <- d3afuns <- rep(NA, nobs)
                d2afuns[keep] <- d2afun(zetas[keep])
                d3afuns[keep] <- d3afun(zetas[keep])
                ## As in brglm2
                s3 <- sum(prior_weights^3 * d3afuns, na.rm = TRUE)
                s2 <- sum(prior_weights^2 * d2afuns, na.rm = TRUE)
                bias_dispersion <- - (nvar - 2) * dispersion^3 / s2  - dispersion^2 * s3 / s2^2
            }
            else {
                bias_dispersion <- NA
            }
        }
        names(bias_dispersion) <- "bias_dispersion"
        names(bias_beta) <- paste0("bias_", names(bias_beta))
        c(bias_beta, bias_dispersion)
    }

    return(list(score = score,
                information = information,
                bias = bias))

}


`compute_auxiliary_functions` <- function(object, ...) {
    UseMethod('compute_auxiliary_functions')
}


`compute_score_mle.glm` <- function(object, ...) {
    object <- enrich(object, with = "auxialiary functions")

}


`compute_score_mle` <- function(object, ...) {
    UseMethod('compute_score_mle')
}


`compute_dispersion_mle.glm` <- function(object, ...) {
    if (object$family$family %in% c("poisson", "binomial")) {
        dispersion_mle <- 1
    }
    else {
        object <- enrich(object, with = "auxiliary functions")
        nobs <- nobs(object)
        prior_weights <- weights(object, type = "prior")
        keep <- prior_weights > 0
        dfResidual <- sum(keep) - object$rank

        gradfun <- function(dispersion) {
            object$auxiliary_functions$score(coef(object), dispersion)[length(coef(object)) + 1]
        }

        if (dfResidual > 0) {
            dispFit <- try(uniroot(f = gradfun, lower = .Machine$double.eps, upper = 10000, tol = 1e-06), silent = FALSE)
            if (inherits(dispFit, "try-error")) {
                warning("the mle of dispersion could not be calculated")
                dispersion_mle <- NA
            }
            else {
                dispersion_mle <- dispFit$root
            }
        }
        else {
            ## if the model is saturated dispersion_mle is NA
            dispersion_mle <- NA
    }
        dispersion_mle
    }
}

`compute_dispersion_mle` <- function(object, ...) {
    UseMethod('compute_dispersion_mle')
}


`compute_expected_information_mle.glm` <- function(object, dispersion = dispersion_mle) {
    object <- enrich(object, with = "auxiliary functions")
    dispersion_mle <- enrich(object, with = "mle of dispersion")$dispersion_mle
    object$auxiliary_functions$information(coef(object), dispersion, type = "expected")
}


`compute_expected_information_mle` <- function(object, ...) {
    UseMethod('compute_expected_information_mle')
}


`compute_observed_information_mle.glm` <- function(object, dispersion = dispersion_mle) {
    object <- enrich(object, with = "auxiliary functions")
    dispersion_mle <- enrich(object, with = "mle of dispersion")$dispersion_mle
    object$auxiliary_functions$information(coef(object), dispersion, type = "observed")
}


`compute_observed_information_mle` <- function(object, ...) {
    UseMethod('compute_observed_information_mle')
}


`compute_bias_mle.glm` <- function(object, ...) {
    ## Write some code to compute the component bias_mle using
    ## the components of object and any other arguments in ...
    ## For example
    cat('component', 'bias_mle', 'for objects of class', 'glm', "\n")
}


`compute_bias_mle` <- function(object, ...) {
    object <- enrich(object, with = "auxiliary functions")
    dispersion_mle <- enrich(object, with = "mle of dispersion")$dispersion_mle
    object$auxiliary_functions$bias(coef(object), dispersion_mle)
}






## ## Call that produced the enrichwith template for the current script:
## create_enrichwith_skeleton(class = "glm", option = c("auxiliary functions",
##     "score vector", "mle of dispersion", "expected information",
##     "observed information", "first-order bias"), description = c("various likelihood-based quantities as functions of the model parameters",
##     "gradient of the log-likelihood at the mle", "mle of the dispersion parameter",
##     "expected information matrix evaluated at the mle", "observed information matrix evaluated at the mle",
##     "first term in the expansion of the bias of the mle at the mle"),
##     component = list("auxiliary_functions", "score_mle", "dispersion_mle",
##         "expected_information_mle", "observed_information_mle",
##         "bias_mle"), path = "~/Downloads", attempt_rename = FALSE)

