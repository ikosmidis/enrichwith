#' Enrich objects of class \code{\link{glm}}
#'
#' Enrich objects of class \code{\link{glm}} with any or all of a set
#' of auxiliary functions, the maximum likelihood estimate of the
#' dispersion parameter, the expected or observed information at the
#' maximum likelihood estimator, and the first term in the expansion
#' of the bias of the maximum likelihood estimator.
#'
#' @param object an object of class glm
#' @param with a character vector of options for the enrichment of \code{object}
#' @param ... extra arguments to be passed to the
#'     \code{compute_*} functions
#'
#' @details
#'
#' The \code{auxiliary_functions} component consists of any or all of the following functions:
#' \itemize{
#' \item \code{score}: the log-likelihood derivatives as a function of the model parameters; see \code{\link{get_score_function.glm}}
#' \item \code{information}: the expected or observed information as a function of the model parameters; see \code{\link{get_information_function.glm}}
#' \item \code{bias}: the first-order term in the expansion of the bias of the maximum likelihood estimator as a function of the model parameters; see \code{\link{get_bias_function.glm}}
#' \item \code{simulate}: a \code{\link{simulate}} function for \code{\link{glm}} objects that can simulate variates from the model at user-supplied parameter values for the regression parameters and the dispersion (default is the maximum likelihood estimates); see \code{\link{get_simulate_function.glm}}
#' \item \code{dmodel}: computes densities or probability mass functions under the model at user-supplied \code{\link{data.frame}}s and at user-supplied values for the regression parameters and the dispersion, if any (default is at the maximum likelihood estimates); see \code{\link{get_dmodel_function.glm}}
#' \item \code{pmodel}: computes distribution functions under the model at user-supplied \code{\link{data.frame}}s and at user-supplied values for the regression parameters and the dispersion, if any (default is at the maximum likelihood estimates); see \code{\link{get_pmodel_function.glm}}
#' \item \code{qmodel}: computes quantile functions under the model at user-supplied \code{\link{data.frame}}s and at user-supplied values for the regression parameters and the dispersion, if any (default is at the maximum likelihood estimates); see \code{\link{get_qmodel_function.glm}}
#' }
#'
#' @return
#'
#' The object \code{object} of class \code{\link{glm}} with extra
#' components. See \code{get_enrichment_options.glm()} for the
#' components and their descriptions.
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
#' clotting <- data.frame(
#'    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
#'    time = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
#'    lot = factor(c(rep(1, 9), rep(2, 9))))
#' cML <- glm(time ~ lot*log(u), data = clotting, family = Gamma)
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
#' biases <- Vectorize(biasfun)(qlogis(probabilities))
#' plot(probabilities, biases, type = "l", ylim = c(-0.5, 0.5),
#'      xlab = expression(pi), ylab = "first-order bias")
#' abline(h = 0, lty = 2); abline(v = 0.5, lty = 2)
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
    out$description <- c('various likelihood-based quantities (gradient of the log-likelihood, expected and observed information matrix and first term in the expansion of the bias of the mle) and a simulate method as functions of the model parameters', 'gradient of the log-likelihood at the mle', 'mle of the dispersion parameter', 'expected information matrix evaluated at the mle', 'observed information matrix evaluated at the mle', 'first term in the expansion of the bias of the mle at the mle')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('auxiliary_functions', 'score_mle', 'dispersion_mle', 'expected_information_mle', 'observed_information_mle', 'bias_mle')
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


`compute_auxiliary_functions.glm` <- function(object, ...) {
    if (is.null(object$model)) {
        object <- update(object, model = TRUE)
    }

    ## Extract formula
    formula <- formula(object)

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
    off <- model.offset(object$model)
    prior_weights <- weights(object, type = "prior")
    keep <- prior_weights > 0
    df_residual <- sum(keep) - object$rank

    ## Take care of aliasing
    na_coefficients <- is.na(coef(object))
    has_na <- any(na_coefficients)

    nvar <- sum(!na_coefficients)

    if (is.null(off)) {
        off <- rep(0, nobs)
    }

    score <- function(coefficients, dispersion, contributions = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (has_na) {
            predictors <- drop(x[, !na_coefficients] %*% coefficients[!na_coefficients] + off)
        }
        else {
            predictors <- drop(x %*% coefficients + off)
        }
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        ## Score for coefficients
        ## score_beta <- colSums(prior_weights * d1mus * (y - fitted_values) * x / variances)/dispersion
        score_beta <- prior_weights * d1mus * (y - fitted_values) * x / variances /dispersion
        ## Score for dispersion
        if (family$family %in% c("poisson", "binomial")) {
            score_dispersion <- NULL
            vnames <- colnames(score_beta)
        }
        else {
            zetas <- -prior_weights/dispersion
            d1afuns <- rep(NA, nobs)
            d1afuns[keep] <- d1afun(zetas[keep])
            if (family$family == "Gamma") d1afuns <- d1afuns - 2
            devianceResiduals <- family$dev.resids(y, fitted_values, prior_weights)
            EdevianceResiduals <- prior_weights * d1afuns
            ## score_dispersion <- sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / (2 * dispersion^2)
            score_dispersion <- (devianceResiduals - EdevianceResiduals) / (2 * dispersion^2)
            vnames <- c(colnames(score_beta), "dispersion")
        }
        if (has_na) {
            score_beta[, na_coefficients] <- NA
        }
        ## Overall score
        out <- cbind(score_beta, score_dispersion)
        colnames(out) <- vnames

        out <- if (contributions) out else colSums(out)
        attr(out, "coefficients") <- coefficients
        attr(out, "dispersion") <- dispersion
        out
    }

    information <- function(coefficients, dispersion,
                            type = c("expected", "observed"), QR = FALSE, CHOL = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (has_na) {
            predictors <- drop(x[, !na_coefficients] %*% coefficients[!na_coefficients] + off)
        }
        else {
            predictors <- drop(x %*% coefficients + off)
        }
        type <- match.arg(type)
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
        if (has_na) {
            info_beta[na_coefficients, ] <- NA
            info_beta[, na_coefficients] <- NA
        }
        ## If there is no dispersion parameter then return the
        ## information for the coefficients only
        if (family$family %in% c("poisson", "binomial")) {
            out <- info_beta
            colnames(out) <- rownames(out) <- colnames(x)
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
                if (family$family == "Gamma") d1afuns <- d1afuns - 2
                devianceResiduals <- family$dev.resids(y, fitted_values, prior_weights)
                EdevianceResiduals <- prior_weights * d1afuns
                info_dispe <- info_dispe + sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE) / dispersion^3
            }
            if (has_na) {
                info_cross[na_coefficients] <- NA
            }
            out <- rbind(cbind(info_beta, info_cross),
                         c(info_cross, info_dispe))
            colnames(out) <- rownames(out) <- c(colnames(x), "dispersion")
        }
        if (CHOL)
            out <- chol(out)
        attr(out, "coefficients") <- coefficients
        attr(out, "dispersion") <- dispersion
        out
    }

    bias <- function(coefficients, dispersion) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (has_na) {
            predictors <- drop(x[, !na_coefficients] %*% coefficients[!na_coefficients] + off)
        }
        else {
            predictors <- drop(x %*% coefficients + off)
        }
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        d2mus <- d2mu(predictors)
        variances <- variance(fitted_values)
        working_weights <- prior_weights * d1mus^2 / variances
        Qr <- information(coefficients, dispersion = dispersion, QR = TRUE)
        inds <- seq.int(Qr$rank)
        Q <- qr.Q(Qr)[, inds, drop = FALSE]
        if (all(dim(Q) == c(1, 1))) {
            hats <- 1
        }
        else {
            hats <- rowSums(Q * Q)
        }
        ksi <- -0.5 * dispersion * d2mus * hats / (d1mus * sqrt(working_weights))

        bias_beta <- numeric(ncol(x))
        ## coefnames <- names(coefficients)
        ## if (is.null(coefnames)) {
        coefnames <- colnames(x)
        ## }
        names(bias_beta) <- coefnames
        biases <- drop(tcrossprod(ksi %*% Q, solve(qr.R(Qr)[inds, inds, drop = FALSE])))
        bias_beta[names(biases)] <- biases
        if (family$family %in% c("poisson", "binomial")) {
            bias_dispersion <- NULL
            vnames <- names(bias_beta)
        }
        else {
            if (df_residual > 0) {
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
        if (has_na) {
            bias_beta[na_coefficients] <- NA
        }
        out <- c(bias_beta, bias_dispersion)
        names(out) <- vnames
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
                               rnorm(nsim * n, mean = fitted_values, sd = sqrt(dispersion/prior_weights))
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
                               else rbinom(n * nsim, size = prior_weights, prob = fitted_values)/prior_weights
                           },
                           "poisson" = {
                               if (any(prior_weights != 1)) {
                                   warning("ignoring prior weights")
                               }
                               rpois(n * nsim, lambda = fitted_values)
                           },
                           "inverse.gaussian" = {
                               SuppDists::rinvGauss(n * nsim, nu = fitted_values, lambda = prior_weights/dispersion)
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

    ## data should have a response
    dmodel <- function(data, coefficients, dispersion, log = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (missing(data)) {
            mf <- object$model
        }
        else  {
            mf <- model.frame(formula = formula, data = data)
        }
        contr <- attr(model.matrix(object), "contrasts")
        new_x <- model.matrix(object = formula, data = mf, terms = terms, contrasts.arg = contr)
        new_y <- model.response(mf)
        new_off <- model.offset(mf)
        if (is.null(new_off)) {
            new_off <- rep(0, nrow(mf))
        }
        if (missing(data)) {
            new_prior_weights <- model.weights(mf)
        }
        else {
            new_prior_weights <- with(data, eval(object$call$weights))
        }
        if (is.null(new_prior_weights)) {
            new_prior_weights <- rep(1, nrow(mf))
        }
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (has_na) {
            predictors <- drop(new_x[, !na_coefficients] %*% coefficients[!na_coefficients] + new_off)
        }
        else {
            predictors <- drop(new_x %*% coefficients + new_off)
        }
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        dfun <- switch(family$family,
                       "gaussian" = {
                           dnorm(new_y, mean = fitted_values, sd = sqrt(dispersion/new_prior_weights), log = log)
                       },
                       "Gamma" = {
                           if (any(new_prior_weights!= 1)) {
                               message("using prior weights in the shape parameters")
                           }
                           dgamma(new_y, shape = new_prior_weights/dispersion, scale = fitted_values*dispersion, log = log)
                       },
                       "binomial" = {
                           if (any(new_prior_weights %% 1 != 0)) {
                               stop("cannot simulate from non-integer prior.weights")
                           }
                           if (is.matrix(new_y) && ncol(new_y)) {
                               new_prior_weights <- rowSums(new_y)
                               new_y <- new_y[, 1]
                               dbinom(new_y, size = new_prior_weights, prob = fitted_values, log = log)
                           }
                           else {
                               if (is.factor(new_y)) {
                                   new_y <- as.numeric(new_y) - 1
                                   dbinom(new_y, size = 1, prob = fitted_values, log = log)
                               }
                               else {
                                   dbinom(new_y*new_prior_weights, size = new_prior_weights, prob = fitted_values, log = log)
                               }
                           }
                       },
                       "poisson" = {
                           if (any(new_prior_weights != 1)) {
                                   warning("ignoring prior weights")
                           }
                           dpois(new_y, lambda = fitted_values, log = log)
                       },
                       "inverse.gaussian" = {
                           SuppDists::dinvGauss(new_y, nu = fitted_values, lambda = new_prior_weights/dispersion, log = log)
                       },
                       NULL)
        attr(dfun, "coefficients") <- coefficients
        attr(dfun, "dispersion") <- dispersion
        dfun
    }

    ## data should have a response
    pmodel <- function(data, coefficients, dispersion, lower.tail = TRUE, log.p = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        contr <- attr(model.matrix(object), "contrasts")
        if (missing(data)) {
            mf <- object$model
        }
        else  {
            mf <- model.frame(formula = formula, data = data)
        }
        new_x <- model.matrix(object = formula, data = mf, terms = terms, contrasts.arg = contr)
        new_y <- model.response(mf)
        new_off <- model.offset(mf)
        if (is.null(new_off)) {
            new_off <- rep(0, nrow(mf))
        }
        if (missing(data)) {
            new_prior_weights <- model.weights(mf)
        }
        else {
            new_prior_weights <- with(data, eval(object$call$weights))
        }
        if (is.null(new_prior_weights)) {
            new_prior_weights <- rep(1, nrow(mf))
        }
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (has_na) {
            predictors <- drop(new_x[, !na_coefficients] %*% coefficients[!na_coefficients] + new_off)
        }
        else {
            predictors <- drop(new_x %*% coefficients + new_off)
        }
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        pfun <- switch(family$family,
                       "gaussian" = {
                           pnorm(new_y, mean = fitted_values, sd = sqrt(dispersion/new_prior_weights), lower.tail = lower.tail, log.p = log.p)
                       },
                       "Gamma" = {
                           if (any(new_prior_weights!= 1)) {
                               message("using prior weights in the shape parameters")
                           }
                           pgamma(new_y, shape = new_prior_weights/dispersion, scale = fitted_values*dispersion, lower.tail = lower.tail, log.p = log.p)
                       },
                       "binomial" = {
                           if (any(new_prior_weights %% 1 != 0)) {
                               stop("cannot simulate from non-integer prior.weights")
                           }
                           if (is.matrix(new_y) && ncol(new_y)) {
                               new_prior_weights <- rowSums(new_y)
                               new_y <- new_y[, 1]
                               pbinom(new_y, size = new_prior_weights, prob = fitted_values, log.p = log.p)
                           }
                           else {
                               if (is.factor(new_y)) {
                                   new_y <- as.numeric(new_y) - 1
                                   pbinom(new_y, size = 1, prob = fitted_values, log.p = log.p)
                               }
                               else {
                                   pbinom(new_y*new_prior_weights, size = new_prior_weights, prob = fitted_values, log.p = log.p)
                               }
                           }
                       },
                       "poisson" = {
                           if (any(new_prior_weights != 1)) {
                                   warning("ignoring prior weights")
                           }
                           ppois(new_y, lambda = fitted_values, log.p = log.p)
                       },
                       "inverse.gaussian" = {
                           SuppDists::pinvGauss(new_y, nu = fitted_values, lambda = new_prior_weights/dispersion, lower.tail = lower.tail, log.p = log.p)
                       },
                       NULL)
        attr(pfun, "coefficients") <- coefficients
        attr(pfun, "dispersion") <- dispersion
        pfun
    }


    ## any response in the data is ignored
    qmodel <- function(p, data, coefficients, dispersion, lower.tail = TRUE, log.p = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        ## output an function that takes as input a data frame and returns densities
        contr <- attr(model.matrix(object), "contrasts")
        if (missing(data)) {
            mf <- object$model
        }
        else  {
            mf <- model.frame(formula = formula, data = data)
        }
        if (length(p) != nrow(mf)) {
            stop("length(p) must be equal to nrow(data)")
        }
        new_x <- model.matrix(object = formula, data = mf, terms = terms, contrasts.arg = contr)
        new_y <- model.response(mf)
        new_off <- model.offset(mf)
        if (is.null(new_off)) {
            new_off <- rep(0, nrow(mf))
        }
        if (missing(data)) {
            new_prior_weights <- model.weights(mf)
        }
        else {
            new_prior_weights <- with(data, eval(object$call$weights))
        }
        if (is.null(new_prior_weights)) {
            new_prior_weights <- rep(1, nrow(mf))
        }
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
        if (has_na) {
            predictors <- drop(new_x[, !na_coefficients] %*% coefficients[!na_coefficients] + new_off)
        }
        else {
            predictors <- drop(new_x %*% coefficients + new_off)
        }
        fitted_values <- linkinv(predictors)
        d1mus <- d1mu(predictors)
        variances <- variance(fitted_values)
        qfun <- switch(family$family,
                       "gaussian" = {
                           qnorm(p, mean = fitted_values, sd = sqrt(dispersion/new_prior_weights), lower.tail = lower.tail, log.p = log.p)
                       },
                       "Gamma" = {c
                           if (any(new_prior_weights!= 1)) {
                               message("using prior weights in the shape parameters")
                           }
                           qgamma(p, shape = new_prior_weights/dispersion, scale = fitted_values*dispersion, lower.tail = lower.tail, log.p = log.p)
                       },
                       "binomial" = {
                           if (any(new_prior_weights %% 1 != 0)) {
                               stop("cannot simulate from non-integer prior.weights")
                           }
                           if (is.matrix(new_y) && ncol(new_y)) {
                               new_prior_weights <- rowSums(new_y)
                               new_y <- new_y[, 1]
                               qbinom(p, size = new_prior_weights, prob = fitted_values, lower.tail = lower.tail, log.p = log.p)
                           }
                           else {
                               if (is.factor(new_y)) {
                                   new_y <- as.numeric(new_y) - 1
                                   qbinom(p, size = 1, prob = fitted_values, lower.tail = lower.tail, log.p = log.p)
                               }
                               else {
                                   qbinom(p, size = new_prior_weights, prob = fitted_values, lower.tail = lower.tail, log.p = log.p)
                               }
                           }
                       },
                       "poisson" = {
                           if (any(new_prior_weights != 1)) {
                                   warning("ignoring prior weights")
                           }
                           qpois(p, lambda = fitted_values, log.p = log.p)
                       },
                       "inverse.gaussian" = {
                           SuppDists::qinvGauss(p, nu = fitted_values, lambda = new_prior_weights/dispersion, lower.tail = lower.tail, log.p = log.p)
                       },
                       NULL)
        attr(qfun, "coefficients") <- coefficients
        attr(qfun, "dispersion") <- dispersion
        qfun
    }


    ## any response in the data is ignored
    ## To be implemented at a later release
    rmodel <- function(n, data, coefficients, dispersion, nsim = 1, seed = NULL) {
        if (missing(coefficients)) {
            coefficients <- coef(object)
        }
        if (missing(dispersion)) {
            dispersion <- enrich(object, with = "mle of dispersion")$dispersion_mle
        }
    }

    return(list(score = score,
                information = information,
                bias = bias,
                simulate = simulate,
                dmodel = dmodel,
                pmodel = pmodel,
                qmodel = qmodel))

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
        df_residual <- sum(keep) - object$rank

        gradfun <- function(logdispersion) {
            beta <- coef(object, model = "mean")
            object$auxiliary_functions$score(beta, exp(logdispersion))[length(beta) + 1]
        }
        if (df_residual > 0) {
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


`compute_expected_information_mle.glm` <- function(object, dispersion) {
    get_information_function(object)(dispersion = dispersion, type = "expected")
}


`compute_expected_information_mle` <- function(object, ...) {
    UseMethod('compute_expected_information_mle')
}


`compute_observed_information_mle.glm` <- function(object, dispersion = dispersion) {
    get_information_function(object)(dispersion = dispersion, type = "observed")
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
               object <- enrich(object, with = "mle of dispersion")
               object$dispersion_mle
           },
           full = {
               object <- enrich(object, with = "mle of dispersion")
               c(beta, object$dispersion_mle)
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
#' \item{QR}{If \code{TRUE}, then the QR decomposition of \deqn{W^{1/2} X} is returned, where \deqn{W} is a diagonal matrix with the working weights (\code{object$weights}) and \deqn{X} is the model matrix.}
#'
#' \item{CHOL}{If \code{TRUE}, then the Cholesky decomposition of the information matrix at the coefficients is returned}
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


#' Function to compute/extract a \code{dmodel} function
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{data}{a data frame with observations at which to compute
#' densities. If missing then densities are computed at the model
#' frame extracted from the object (see \code{\link{glm}})}
#'
#' \item{coefficients}{the regression coefficients at which the
#' densities are computed. If missing then the maximum likelihood
#' estimates are used}
#'
#' \item{dispersion}{the dispersion parameter at which the densities
#' function is computed. If missing then the maximum likelihood
#' estimate is used}
#'
#' \item{log}{logical; if \code{TRUE}, the logarithm of the density is
#' returned}
#'
#' }
#'
#' @export
get_dmodel_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$dmodel
    }
    else {
        object$auxiliary_functions$dmodel
    }
}

#' Function to compute/extract a \code{pmodel} function
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{data}{a data frame with observations at which to compute the
#' distribution function. If missing then probabilities are computed
#' at the model frame extracted from the object (see
#' \code{\link{glm}})}
#'
#' \item{coefficients}{the regression coefficients at which the
#' distribution function are computed. If missing then the maximum
#' likelihood estimates are used}
#'
#' \item{dispersion}{the dispersion parameter at which the
#' distribution function is computed. If missing then the maximum
#' likelihood estimate is used}
#'
#' \item{log.p}{logical; if \code{TRUE}, the logarithm of the
#' distribution function is returned}
#'
#' \item{lower.tail}{logical; if \code{TRUE} (default), probabilities
#' are P[X <= x] otherwise, P[X > x]}
#' }
#'
#' @export
get_pmodel_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$pmodel
    }
    else {
        object$auxiliary_functions$pmodel
    }
}

#' Function to compute/extract a \code{qmodel} function
#'
#' @param object an object of class \code{glm} or\code{enriched_glm}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{p}{a vector of probabilities with \code{length(p)} equal to
#' \code{nrow(data)} at which to evaluate quantiles}
#'
#' \item{data}{a data frame with observations at which to compute the
#' quantiles. If missing then quantiles are computed at the model
#' frame extracted from the object (see \code{\link{glm}})}
#'
#' \item{coefficients}{the regression coefficients at which the
#' quantiles are computed. If missing then the maximum likelihood
#' estimates are used}
#'
#' \item{dispersion}{the dispersion parameter at which the
#' quantiles are computed. If missing then the maximum
#' likelihood estimate is used}
#'
#' \item{log.p}{logical; if \code{TRUE}, the logarithm of the
#' probabilities is used}
#'
#' \item{lower.tail}{logical; if \code{TRUE} (default), probabilities
#' are P[X <= x] otherwise, P[X > x]}
#'
#' }
#'
#' @export
get_qmodel_function.glm <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$qmodel
    }
    else {
        object$auxiliary_functions$qmodel
    }
}


## ## ## Call that produced the enrichwith template for the current script:
## create_enrichwith_skeleton(class = "glm", option = c("auxiliary functions",
##     "score vector", "mle of dispersion", "expected information",
##     "observed information", "first-order bias"), description = c("various likelihood-based quantities (gradient of the log-likelihood, expected and observed information matrix and first term in the expansion of the bias of the mle) and a simulate method as functions of the model parameters",
##     "gradient of the log-likelihood at the mle", "mle of the dispersion parameter",
##     "expected information matrix evaluated at the mle", "observed information matrix evaluated at the mle",
##     "first term in the expansion of the bias of the mle at the mle"),
##     component = list("auxiliary_functions", "score_mle", "dispersion_mle",
##         "expected_information_mle", "observed_information_mle",
##         "bias_mle"), path = "~/Downloads", attempt_rename = FALSE)
