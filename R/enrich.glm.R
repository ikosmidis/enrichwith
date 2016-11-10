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
#' expected or observed information, the first-order bias of the
#' maximum likelihood estimator as functions of the model parameters,
#' and a \code{simulate} function that takes as input the model
#' parameters (including the dispersion if any). The result from the
#' \code{simulate} auxiliary function has the same structure to that
#' of the \code{\link{simulate}} method for \code{\link{glm}} objects.
#'
#' @return
#'
#' The object \code{object} of class \code{\link{glm}} with extra
#' components. \code{get_enrichment_options.glm()} returns the
#' components and their descriptions.
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
#' clotting <- data.frame(
#'    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
#'    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
#'    lot = factor(c(rep(1, 9), rep(2, 9))))
#' cML <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)
#'
#' # The simulate method for the above fit would simulate at coef(cML)
#' # for the regression parameters and MASS::gamma.dispersion(cML) for
#' # the dispersion. It is not possible to simulate at different
#' # parameter values than those, at least not, without "hacking" the
#' # cML object.
#'
#' # A general simulator for cML results via its enrichment with
#' # auxiliary functions:
#' cML_functions <- get_auxiliary_functions(cML)
#' # which is a shorthand for
#' # enriched_cML <- enrich(cML, with = "auxiliary functions")
#' # cML_functions <- enriched_cML$auxiliary_functions
#'
#' # To simulate 2 samples at the maximum likelihood estimator do
#' dispersion_mle <- MASS::gamma.dispersion(cML)
#' cML_functions$simulate(coef = coef(cML),
#'                        dispersion = dispersion_mle,
#'                        nsim = 2, seed = 123)
#' # To simulate 5 samples at c(0.1, 0.1, 0, 0) and dispersion 0.2 do
#' cML_functions$simulate(coef = c(0.1, 0.1, 0, 0),
#'                        dispersion = 0.2,
#'                        nsim = 5, seed = 123)
#'
#' }
#'
#' \dontrun{
#'
#' ## Reproduce left plot in Figure 4.1 in Kosimdis (2007)
#' ## (see http://www.ucl.ac.uk/~ucakiko/files/ikosmidis_thesis.pdf)
#' mod <- glm(1 ~ 1, weights = 10, family = binomial())
#' enriched_mod <- enrich(mod, with = "auxiliary functions")
#' biasfun <- enriched_mod$auxiliary_functions$bias
#' probabilities <- seq(1e-02, 1 - 1e-02, length = 100)
#' biases <- Vectorize(biasfun)(qlogis(probabilities))[1,]
#' plot(probabilities, biases, type = "l", ylim = c(-0.5, 0.5),
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
    out$description <- c('various likelihood-based quantities (gradient of the log-likelihood, expected and observed information matrix and first term in the expansion of the bias of the mle) as functions of the model parameters', 'gradient of the log-likelihood at the mle', 'mle of the dispersion parameter', 'expected information matrix evaluated at the mle', 'observed information matrix evaluated at the mle', 'first term in the expansion of the bias of the mle at the mle')
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
    y <- object$y
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

    score <- function(coefficients, dispersion) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        predictors <- drop(x %*% coefficients + off)
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        ## Score for coefficients
        score_beta <- colSums(prior_weights * d1mus * (y - fitted_values) * x / variances)/dispersion
        ## Score for dispersion
        if (family$family %in% c("poisson", "binomial")) {
            score_dispersion <- NULL
            vnames <- names(score_beta)
        }
        else {
            zetas <- -prior_weights/dispersion
            d1afuns <- rep(NA, nobs)
            d1afuns[keep] <- d1afun(zetas[keep])
            devianceResiduals <- family$dev.resids(y, fitted_values, prior_weights)
            EdevianceResiduals <- prior_weights * d1afuns
            score_dispersion <- sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / (2 * dispersion^2)
            vnames <- c(names(score_beta), "dispersion")
        }
        ## Overall score
        out <- c(score_beta, score_dispersion)
        names(out) <- paste0("grad_", vnames)
        attr(out, "coefficients") <- coefficients
        attr(out, "dispersion") <- dispersion
        out
    }

    information <- function(coefficients, dispersion,
                            type = c("expected", "observed"), QR = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        type <- match.arg(type)
        predictors <- drop(x %*% coefficients + off)
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        working_weights <- prior_weights * d1mus^2 / variances
        wx <- x * sqrt(working_weights)
        if (QR) {
            return(qr(wx))
        }
        ## expected info coefficients-coefficients
        info_beta <- crossprod(wx) / dispersion
        if (type == "observed") {
            d2mus <- d2mu(predictors)
            d1variances <- d1variance(fitted_values)
            w <- prior_weights * (d2mus / variances - d1mus^2 * d1variances / variances^2) * (y - fitted_values)
            ## observed info coefficients-coefficients
            info_beta <- info_beta - t(x * w) %*% x / dispersion
        }
        rownames(info_beta) <- colnames(info_beta) <- colnames(x)
        ## If there is no dispersion parameter then return the
        ## information for the coefficients only
        if (family$family %in% c("poisson", "binomial")) {
            return(info_beta)
        }
        ## If there is a dispersion parameter then return the
        ## information on the coefficients and the dispersion
        else {
            ## expected info coefficients-dispersion
            info_cross <- rep(0, ncol(info_beta))
            ## expected info dispersion-dispersion
            zetas <- -prior_weights/dispersion
            d2afuns <- rep(NA, nobs)
            d2afuns[keep] <- d2afun(zetas[keep])
            info_dispe <- sum(prior_weights^2 * d2afuns, na.rm = TRUE) / (2 * dispersion^4)
            if (type == "observed") {
                ## observed info coefficients-dispersion
                info_cross <- info_cross + colSums(prior_weights * d1mus * (y - fitted_values) * x / variances)/dispersion^2
                ## observed info dispersion-dispersion
                zetas <- -prior_weights/dispersion
                d1afuns <- rep(NA, nobs)
                d1afuns[keep] <- d1afun(zetas[keep])
                devianceResiduals <- family$dev.resids(y, fitted_values, prior_weights)
                EdevianceResiduals <- prior_weights * d1afuns
                info_dispe <- info_dispe + sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / dispersion^3
            }
            out <- rbind(cbind(info_beta, info_cross),
                         c(info_cross, info_dispe))
            colnames(out) <- rownames(out) <- c(colnames(x), "dispersion")
            attr(out, "coefficients") <- coefficients
            attr(out, "dispersion") <- dispersion
            out
        }
    }

    bias <- function(coefficients, dispersion) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        predictors <- drop(x %*% coefficients + off)
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        d2mus <- d2mu(predictors)
        variances <- variance(fitted_values)
        working_weights <- prior_weights * d1mus^2 / variances
        Qr <- information(coefficients, dispersion = dispersion, QR = TRUE)
        Q <- qr.Q(Qr)
        hats <- rowSums(Q * Q)
        ksi <- -0.5 * dispersion * d2mus * hats / (d1mus * sqrt(working_weights))
        bias_beta <- drop(tcrossprod(ksi %*% Q, solve(qr.R(Qr))))
        if (family$family %in% c("poisson", "binomial")) {
            bias_dispersion <- NULL
            vnames <- names(bias_beta)
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
            vnames <- c(names(bias_beta), "dispersion")
        }
        out <- c(bias_beta, bias_dispersion)
        names(out) <- paste0("bias_", vnames)
        attr(out, "coefficients") <- coefficients
        attr(out, "dispersion") <- dispersion
        out
    }

    simulate <- function(coefficients, dispersion, nsim = 1, seed = NULL) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        if (is.null(seed))
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        predictors <- drop(x %*% coefficients + off)
        fitted_values <- linkinv(predictors)
        fitted_names <- names(fitted_values)
        n <- length(fitted_values)
        variates <- switch(family$family,
                           "gaussian" = {
                               rnorm(nsim * n, mean = fitted_values, sd = dispersion/prior_weights)
                           },
                           "Gamma" = {
                               if (any(prior_weights!= 1)) {
                                   message("using prior weights in the shape parameters")
                               }
                               rgamma(nsim * n, shape = prior_weights/dispersion, scale = fitted_values*dispersion)
                           },
                           "binomial" = {
                               if (any(prior_weights %% 1 != 0))
                                   stop("cannot simulate from non-integer prior.weights")
                               if (!is.null(mf <- object$model)) {
                                   y <- model.response(mf)
                                   if (is.factor(y)) {
                                       yy <- factor(1 + rbinom(n * nsim, size = 1, prob = fitted_values),
                                                    labels = levels(y))
                                       split(yy, rep(seq_len(nsim), each = n))
                                   }
                                   else if (is.matrix(y) && ncol(y) == 2) {
                                       yy <- vector("list", nsim)
                                       for (i in seq_len(nsim)) {
                                           Y <- rbinom(n, size = prior_weights, prob = fitted_values)
                                           YY <- cbind(Y, prior_weights - Y)
                                           colnames(YY) <- colnames(y)
                                           yy[[i]] <- YY
                                       }
                                       yy
                                   }
                                   else rbinom(n * nsim, size = prior_weights, prob = fitted_values)/prior_weights
                               }
                               else rbinom(n & nsim, size = prior_weights, prob = fitted_values)/prior_weights
                           },
                           "poisson" = {
                               if (any(prior_weights != 1)) {
                                   warning("ignoring prior weights")
                               }
                               rpois(nsim * n, lambda = fitted_values)
                           },
                           "inverse.gaussian" = {
                               SuppDists::rinvGauss(nsim * n, nu = fitted_values, lambda = prior_weights/dispersion)
                           },
NULL)
        ## Inspired by stats:::simulate.lm
        if (!is.list(variates)) {
            dim(variates) <- c(n, nsim)
            variates <- as.data.frame(variates)
        }
        else {
            class(variates) <- "data.frame"
        }
        names(variates) <-  paste("sim", seq_len(nsim), sep = "_")
        if (!is.null(fitted_names)) {
            row.names(variates) <- fitted_names
        }
        attr(variates, "seed") <- RNGstate
        attr(variates, "coefficients") <- coefficients
        attr(variates, "dispersion") <- dispersion
        variates
    }

    return(list(score = score,
                information = information,
                bias = bias,
                simulate = simulate))

}


`compute_auxiliary_functions` <- function(object, ...) {
    UseMethod('compute_auxiliary_functions')
}


`compute_score_mle.glm` <- function(object, ...) {
    get_score_function(object)()
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

        gradfun <- function(logdispersion) {
            beta <- coef(object, model = "mean")
            object$auxiliary_functions$score(beta, exp(logdispersion))[length(beta) + 1]
        }

        if (dfResidual > 0) {
            dispFit <- try(uniroot(f = gradfun, lower = 0.5*log(.Machine$double.eps), upper = 20, tol = 1e-08, maxiter = 10000), silent = FALSE)
            if (inherits(dispFit, "try-error")) {
                warning("the mle of dispersion could not be calculated")
                dispersion_mle <- NA
            }
            else {
                dispersion_mle <- exp(dispFit$root)
            }
        }
        else {
            ## if the model is saturated dispersion_mle is NA
            dispersion_mle <- NA
    }
        names(dispersion_mle) <- "dispersion"
        dispersion_mle
    }
}

`compute_dispersion_mle` <- function(object, ...) {
    UseMethod('compute_dispersion_mle')
}


`compute_expected_information_mle.glm` <- function(object, dispersion = dispersion_mle) {
    get_information_function(object)(type = "expected")
}


`compute_expected_information_mle` <- function(object, ...) {
    UseMethod('compute_expected_information_mle')
}


`compute_observed_information_mle.glm` <- function(object, dispersion = dispersion_mle) {
    get_information_function(object)(type = "observed")
}


`compute_observed_information_mle` <- function(object, ...) {
    UseMethod('compute_observed_information_mle')
}


`compute_bias_mle.glm` <- function(object, ...) {
    get_bias_function(object)()
}


`compute_bias_mle` <- function(object, ...) {
    UseMethod('compute_bias_mle')
}


## Other extractor functions

#' Function to extract model coefficients from objects of class \code{enriched_glm}
#'
#' @param object an object of class \code{enriched_glm}
#' @param model either "mean" for the estimates of the parameters in the linear predictor, or "dispersion" for the estimate of the dispersion, or "full" for all estimates
#' @param ... currently unused
#' @export
coef.enriched_glm <- function(object, model = c("mean", "full", "dispersion"), ...) {
    beta <- object$coefficients
    switch(match.arg(model),
           mean = {
               beta
           },
           dispersion = {
               object$dispersion
           },
           full = {
               c(beta, object$dispersion)
           })
}

#' Function to compute/extract auxiliary functions from objects of
#' class \code{glm}/\code{enriched_glm}
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#'
#' See \code{\link{enrich.glm}} for details.
#'
#' @export
get_auxiliary_functions.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        enriched_object <- enrich(object, with = "auxiliary functions")
        enriched_object$auxiliary_functions
    }
    else {
        object$auxiliary_functions
    }
}

#' Function to compute/extract a simulate function for response
#' vectors from an object of class \code{glm}/\code{enriched_glm}
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted simulate function has arguments
#' \describe{
#'
#' \item{coefficients}{the regression coefficients at which the
#' response vectors are simulated. If missing then the maximum
#' likelihood estimates are used}
#'
#' \item{dispersion}{the dispersion parameter at which the response
#' vectors are simulated. If missing then the maximum likelihood
#' estimate is used}
#'
#' \item{nsim}{number of response vectors to simulate.  Defaults to \code{1}}
#'
#' \item{seed}{an object specifying if and how the random number
#' generator should be initialized ('seeded'). It can be either
#' \code{NULL} or an integer that will be used in a call to
#' \code{set.seed} before simulating the response vectors.  If set,
#' the value is saved as the \code{seed} attribute of the returned
#' value.  The default, \code{NULL} will not change the random
#' generator state, and return \code{.Random.seed} as the \code{seed}
#' attribute, see \code{Value}}
#'
#' }
#'
#' @export
get_simulate_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$simulate
    }
    else {
        object$auxiliary_functions$simulate
    }
}


#' Function to compute/extract a function that returns the scores
#' (derivatives of the log-likelihood) for an object of class
#' \code{glm}/\code{enriched_glm}
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{coefficients}{the regression coefficients at which the scores
#' are computed. If missing then the maximum likelihood estimates are
#' used}
#'
#' \item{dispersion}{the dispersion parameter at which the score
#' function is evaluated. If missing then the maximum likelihood
#' estimate is used}
#'
#' }
#'
#' @export
get_score_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$score
    }
    else {
        object$auxiliary_functions$score
    }
}

#' Function to compute/extract a function that returns the information
#' matrix for an object of class \code{glm}/\code{enriched_glm}
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{coefficients}{the regression coefficients at which the
#' information matrix is evaluated. If missing then the maximum
#' likelihood estimates are used}
#'
#' \item{dispersion}{the dispersion parameter at which the information
#' matrix is evaluated. If missing then the maximum likelihood estimate
#' is used}
#'
#' \item{type}{should the function return th 'expected' or 'observed' information? Default is \code{expected}}
#'
#' \item{QR}{If \code{TRUE}, then the QR decomposition of the expected information for the coefficients is returned}
#'
#' }
#'
#' @export
get_information_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$information
    }
    else {
        object$auxiliary_functions$information
    }
}


#' Function to compute/extract a function that returns the first term
#' in the expansion of the bias of the MLE for the parameters of an
#' object of class \code{glm}/\code{enriched_glm}
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{coefficients}{the regression coefficients at which the
#' first-order bias is evacuated. If missing then the maximum
#' likelihood estimates are used}
#'
#' \item{dispersion}{the dispersion parameter at which the first-order
#' bias is evaluated. If missing then the maximum likelihood estimate
#' is used}
#'
#' }
#'
#' @export
get_bias_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$bias
    }
    else {
        object$auxiliary_functions$bias
    }
}




## ## Call that produced the enrichwith template for the current script:
## create_enrichwith_skeleton(class = "glm", option = c("auxiliary functions",
##     "score vector", "mle of dispersion", "expected information",
##     "observed information", "first-order bias"), description = c("various likelihood-based quantities (gradient of the log-likelihood, expected and observed information matrix and first term in the expansion of the bias of the mle) as functions of the model parameters",
##     "gradient of the log-likelihood at the mle", "mle of the dispersion parameter",
##     "expected information matrix evaluated at the mle", "observed information matrix evaluated at the mle",
##     "first term in the expansion of the bias of the mle at the mle"),
##     component = list("auxiliary_functions", "score_mle", "dispersion_mle",
##         "expected_information_mle", "observed_information_mle",
##         "bias_mle"), path = "~/Downloads", attempt_rename = FALSE)
