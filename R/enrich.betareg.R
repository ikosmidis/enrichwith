#' Enrich objects of class betareg
#'
#' Enrich objects of class \code{\link[betareg]{betareg}} with any or all of a
#' set of auxiliary functions, the expected information at the maximum
#' likelihood estimator, and the first term in the expansion of the
#' bias of the maximum likelihood estimator.
#'
#'
#' @param object an object of class \code{\link[betareg]{betareg}}
#' @param with  a character vector of options for the enrichment of \code{object}
#' @param ... extra arguments to be passed to the
#'     \code{compute_*} functions
#'
#' @details
#' The \code{auxiliary_functions} component consists of any or all of the following functions:
#' \itemize{
#' \item \code{score}: the log-likelihood derivatives as a function of the model parameters; see \code{get_score_function.betareg}
#' \item \code{information}: the expected information as a function of the model parameters; see \code{\link{get_information_function.betareg}}
#' \item \code{bias}: the first-order term in the expansion of the bias of the maximum likelihood estimator as a function of the model parameters; see \code{\link{get_bias_function.betareg}}
#' \item \code{simulate}: a \code{\link{simulate}} function for \code{\link[betareg]{betareg}} objects that can simulate variates from the model at user-supplied parameter values for the regression parameters (default is the maximum likelihood estimates); see \code{\link{get_simulate_function.betareg}}
#' }
#'
#' @return The object \code{object} of class \code{\link[betareg]{betareg}}
#'     with extra components. \code{get_enrichment_options.betareg()}
#'     returns the components and their descriptions.
#'
#' @export
#' @examples
#'
#' \dontrun{
#' if (require("betareg")) {
#'
#'    data("GasolineYield", package = "betareg")
#'    gy <- betareg(yield ~ batch + temp, data = GasolineYield)
#'
#'    # Get a function that evaluates the expected information for gy at supplied parameter values
#'    gy_info <- get_information_function(gy)
#'.   # compare standard errors with what `summary` returns
#'    all.equal(sqrt(diag(solve(gy_info())))[1:11],
#'              coef(summary(gy))$mean[, 2], check.attributes = FALSE)
#'.   # evaluating at different parameter values
#'    gy_info(rep(1, length = 12))
#'
#'    # Get a function that evaluates the first-order bias of gy at supplied parameter values
#'    gy_bias <- get_bias_function(gy)
#'    # compare with internal betareg implementation of bias correction
#'    gy_bc <- update(gy, type = "BC")
#'    all.equal(gy_bias(coef(gy)), gy_bc$bias, check.attributes = FALSE)
#'
#'  }
#'}
#'
`enrich.betareg` <- function(object, with = "all", ...) {
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

#' Available options for the enrichment objects of class betareg
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
#' get_enrichment_options.betareg(option = "all")
#' get_enrichment_options.betareg(all_options = TRUE)
#' }
#' @export
`get_enrichment_options.betareg` <- function(object, option, all_options = missing(option)) {
    ## List the enrichment options that you would like to make
    ## available for objects of class
    out <- list()
    out$option <- c('auxiliary functions', 'score vector', 'expected information', 'first-order bias')
    ## Provide the descriptions of the enrichment options
    out$description <- c('various likelihood-based quantities (gradient of the log-likelihood, expected information matrix and first term in the expansion of the bias of the mle) and a simulate method as functions of the model parameters', 'gradient of the log-likelihood at the mle', 'expected information matrix evaluated at the mle', 'first term in the expansion of the bias of the mle at the mle')
    ## Add all as an option
    out$option <- c(out$option, 'all')
    out$description <- c(out$description, 'all available options')
    out$component <- list('auxiliary_functions', 'score_mle', 'expected_information_mle', 'bias_mle')
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


`compute_auxiliary_functions.betareg` <- function(object, ...) {
    if (is.null(object$model)) {
        object <- update(object, model = TRUE)
    }
    y <-  object$y
    x <- model.matrix(object, model = "mean")
    z <- model.matrix(object, model = "precision")
    n <- NROW(x)
    k <- NCOL(x)
    weights <- object$weights
    if (is.null(weights))
        weights <- rep.int(1, n)
    nobs <- sum(weights > 0)
    offset <- object$offset
    offset <- lapply(offset, function(off) {
        if (is.null(off)) rep.int(0, n) else off
    })
    m <- NCOL(z)
    if (m < 1L)
        stop("dispersion regression needs to have at least one parameter")
    phi_const <- (m == 1L) && isTRUE(all.equal(as.vector(z[, 1L]), rep.int(1, n)))
    linkmean <- enrich(object$link$mean)
    linkprec <- enrich(object$link$precision)
    linkinv <- linkmean$linkinv
    mu.eta <- linkmean$mu.eta
    dmu.deta <- linkmean$dmu.deta
    phi_linkinv <- linkprec$linkinv
    phi_mu.eta <- linkprec$mu.eta
    phi_dmu.deta <- linkprec$dmu.deta
    ystar <- qlogis(y)
    u <- log(1 - y)
    score <- function(coefficients, contributions = FALSE) {
        if (missing(coefficients)) {
            coefficients <- coef(object, model = "full")
        }
        beta <- coefficients[seq.int(length.out = k)]
        gamma <- coefficients[seq.int(length.out = m) + k]
        eta <- as.vector(x %*% beta + offset[[1L]])
        phi_eta <- as.vector(z %*% gamma + offset[[2L]])
        mu <- linkinv(eta)
        phi <- phi_linkinv(phi_eta)
        mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
        psi1 <- trigamma(mu * phi)
        psi2 <- trigamma((1 - mu) * phi)
        tbar_ubar <- ystar - mustar
        ubar <- u - digamma((1 - mu) * phi) + digamma(phi)
        rval <- cbind(phi * tbar_ubar * mu.eta(eta) * weights * x,
        (mu * tbar_ubar + ubar) * phi_mu.eta(phi_eta) * weights * z)
        if (contributions)
            rval
        else
            colSums(rval)
    }

    information <- function(coefficients, QR = TRUE, CHOL = FALSE,
                            type = c("expected", "observed")) {
        if (missing(coefficients)) {
            coefficients <- coef(object, model = "full")
        }
        type <- match.arg(type)
        beta <- coefficients[seq.int(length.out = k)]
        gamma <- coefficients[seq.int(length.out = m) + k]
        eta <- as.vector(x %*% beta + offset[[1L]])
        phi_eta <- as.vector(z %*% gamma + offset[[2L]])
        mu <- linkinv(eta)
        phi <- phi_linkinv(phi_eta)
        mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
        psi1 <- trigamma(mu * phi)
        psi2 <- trigamma((1 - mu) * phi)
        a <- psi1 + psi2
        b <- psi1 * mu^2 + psi2 * (1 - mu)^2 - trigamma(phi)
        D1 <- mu.eta(eta)
        D2 <- phi_mu.eta(phi_eta)
        wbb <- phi^2 * a * D1^2
        wpp <- b * D2^2
        wbp <- phi * (mu * a - psi2) * D1 * D2
        kbb <- if (k > 0L)
            crossprod(sqrt(weights) * sqrt(wbb) * x)
        else crossprod(x)
        kpp <- if (m > 0L)
            crossprod(sqrt(weights) * sqrt(wpp) * z)
        else crossprod(z)
        kbp <- if (k > 0L & m > 0L)
            crossprod(weights * wbp * x, z)
        else crossprod(x, z)
        out <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
        if (type == "observed") {
            tbar_ubar <- ystar - mustar
            ubar <- u - digamma((1 - mu) * phi) + digamma(phi)
            D1dash <- dmu.deta(eta)
            D2dash <- phi_dmu.deta(phi_eta)
            kbb <- if (k > 0L)
                       crossprod(phi * D1dash * tbar_ubar * x, x)
                   else crossprod(x)
            kpp <- if (m > 0L)
                       crossprod(D2dash * (mu * tbar_ubar + ubar) * z, z)
                   else crossprod(z)
            kbp <- if (k > 0L & m > 0L)
                       crossprod(D1 * D2 * tbar_ubar * x, z)
                   else crossprod(x, z)
            out <- out - cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
        }
        colnames(out) <- names(coef(object))
        if (CHOL)
            out <- chol(out)
        attr(out, "coefficients") <- coefficients
        out

    }

    bias <- function(coefficients) {
        if (missing(coefficients)) {
            coefficients <- coef(object, model = "full")
        }
        beta <- coefficients[seq.int(length.out = k)]
        gamma <- coefficients[seq.int(length.out = m) + k]
        eta <- as.vector(x %*% beta + offset[[1L]])
        phi_eta <- as.vector(z %*% gamma + offset[[2L]])
        mu <- linkinv(eta)
        phi <- phi_linkinv(phi_eta)
        mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
        psi1 <- trigamma(mu * phi)
        psi2 <- trigamma((1 - mu) * phi)
        InfoInv <- try(solve(information(coefficients)), silent = TRUE)
        D1 <- mu.eta(eta)
        D2 <- phi_mu.eta(phi_eta)
        D1dash <- dmu.deta(eta)
        D2dash <- phi_dmu.deta(phi_eta)
        dPsi1 <- psigamma(mu * phi, 2)
        dPsi2 <- psigamma((1 - mu) * phi, 2)
        kappa2 <- psi1 + psi2
        kappa3 <- dPsi1 - dPsi2
        psi3 <- psigamma(phi, 1)
        dPsi3 <- psigamma(phi, 2)
        PQsum <- function(t) {
            if (t <= k) {
                Xt <- x[, t]
                bb <- if (k > 0L)
                  crossprod(x, weights * phi^2 * D1 * (phi *
                    D1^2 * kappa3 + D1dash * kappa2) * Xt * x)
                else crossprod(x)
                bg <- if ((k > 0L) & (m > 0L))
                  crossprod(x, weights * phi * D1^2 * D2 * (mu *
                    phi * kappa3 + phi * dPsi2 + kappa2) * Xt *
                    z)
                else crossprod(x, z)
                gg <- if (m > 0L)
                  crossprod(z, weights * phi * D1 * D2^2 * (mu^2 *
                    kappa3 - dPsi2 + 2 * mu * dPsi2) * Xt * z) +
                    crossprod(z, weights * phi * D1 * D2dash *
                      (mu * kappa2 - psi2) * Xt * z)
                else crossprod(z)
            }
            else {
                Zt <- z[, t - k]
                bb <- if (k > 0L)
                  crossprod(x, weights * phi * D2 * (phi * D1^2 *
                    mu * kappa3 + phi * D1^2 * dPsi2 + D1dash *
                    mu * kappa2 - D1dash * psi2) * Zt * x)
                else crossprod(x)
                bg <- if ((k > 0L) & (m > 0L))
                  crossprod(x, weights * D1 * D2^2 * (phi * mu^2 *
                    kappa3 + phi * (2 * mu - 1) * dPsi2 + mu *
                    kappa2 - psi2) * Zt * z)
                else crossprod(x, z)
                gg <- if (m > 0L)
                  crossprod(z, weights * D2^3 * (mu^3 * kappa3 +
                    (3 * mu^2 - 3 * mu + 1) * dPsi2 - dPsi3) *
                    Zt * z) + crossprod(z, weights * D2dash *
                    D2 * (mu^2 * kappa2 + (1 - 2 * mu) * psi2 -
                    psi3) * Zt * z)
                else crossprod(z)
            }
            pq <- rbind(cbind(bb, bg), cbind(t(bg), gg))
            sum(diag(InfoInv %*% pq))/2
        }
        if (inherits(InfoInv, "try-error")) {
            bias <- rep.int(NA_real_, k + m)
        }
        else {
            bias <- drop(-InfoInv %*% sapply(1:(k + m), PQsum))
        }
        bias
    }

    simulate <- function(coefficients, nsim = 1, seed = NULL) {
        if (missing(coefficients)) {
            coefficients <-  coef(object, model = "full")
        }
        else {
            if (!isTRUE(identical(length(coefficients), length(coef(object, model = "full"))))) {
                stop("`coefficients` does not have the right length")
            }
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
        beta <- coefficients[seq.int(length.out = k)]
        gamma <- coefficients[seq.int(length.out = m) + k]
        eta <- as.vector(x %*% beta + offset[[1L]])
        phi_eta <- as.vector(z %*% gamma + offset[[2L]])
        mu <- linkinv(eta)
        phi <- phi_linkinv(phi_eta)
        ps <- mu * phi
        qs <- (1 - mu) * phi
        matrix(1 - rbeta(n * nsim, qs, ps), n, nsim)
    }
    return(list(score = score,
                information = information,
                bias = bias,
                simulate = simulate))
}


`compute_auxiliary_functions` <- function(object, ...) {
    UseMethod('compute_auxiliary_functions')
}


`compute_score_mle.betareg` <- function(object, ...) {
    get_score_function(object)()
}


`compute_score_mle` <- function(object, ...) {
    UseMethod('compute_score_mle')
}

`compute_expected_information_mle.betareg` <- function(object, ...) {
    get_information_function(object)()
}


`compute_expected_information_mle` <- function(object, ...) {
    UseMethod('compute_expected_information_mle')
}

`compute_bias_mle.betareg` <- function(object, ...) {
    get_bias_function(object)()
}


`compute_bias_mle` <- function(object, ...) {
    UseMethod('compute_bias_mle')
}

#' Function to compute/extract auxiliary functions from objects of
#' class \code{betreg}/\code{enriched_betareg}
#'
#' @param object an object of class \code{betareg} or\code{enriched_betareg}
#' @param ... currently not used
#'
#' @details
#'
#' See \code{\link{enrich.betareg}} for details.
#'
#' @export
get_auxiliary_functions.betareg <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        enriched_object <- enrich(object, with = "auxiliary functions")
        enriched_object$auxiliary_functions
    }
    else {
        object$auxiliary_functions
    }
}

#' Function to compute/extract a simulate function for response
#' vectors from an object of class \code{betareg}/\code{enriched_betareg}
#'
#' @param object an object of class \code{betareg} or\code{enriched_betareg}
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
get_simulate_function.betareg <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$simulate
    }
    else {
        object$auxiliary_functions$simulate
    }
}

#' Function to compute/extract a function that returns the scores
#' (derivatives of the log-likelihood) for an object of class
#' \code{betareg}/\code{enriched_betareg}
#'
#' @param object an object of class \code{betareg} or\code{enriched_betareg}
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
#' }
#'
#' @export
get_score_function.betareg <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$score
    }
    else {
        object$auxiliary_functions$score
    }
}

#' Function to compute/extract a function that returns the information
#' matrix for an object of class \code{betareg}/\code{enriched_betareg}
#'
#' @param object an object of class \code{betareg} or\code{enriched_betareg}
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
#'
#' \item{type}{should the function return th 'expected' or 'observed' information? Default is \code{expected}}
#'
#' \item{QR}{Currently not used}
#'
#' \item{CHOL}{If \code{TRUE}, then the Cholesky decomposition of the information matrix at the coefficients is returned}
#'
#' }
#'
#' @export
get_information_function.betareg <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$information
    }
    else {
        object$auxiliary_functions$information
    }
}


#' Function to compute/extract a function that returns the first term
#' in the expansion of the bias of the MLE for the parameters of an
#' object of class \code{betareg}/\code{enriched_betareg}
#'
#' @param object an object of class \code{betareg} or\code{enriched_betareg}
#' @param ... currently not used
#'
#' @details
#' The computed/extracted function has arguments
#' \describe{
#'
#' \item{coefficients}{the regression coefficients at which the
#' first-order bias is evacuated. If missing then the maximum
#' likelihood estimates are used}
#' }
#'
#' @export
get_bias_function.betareg <- function(object, ...) {
    if (is.null(object$auxiliary_functions)) {
        get_auxiliary_functions(object)$bias
    }
    else {
        object$auxiliary_functions$bias
    }
}



## ## Call that produced the original version of the enrichwith template for the current script:
## enrichwith:::create_enrichwith_skeleton(class = "betareg", option = c("auxiliary functions",
##     "score vector", "mle of dispersion", "expected information",
##     "first-order bias"), description = c("various likelihood-based quantities (gradient of the log-likelihood, expected information matrix and first term in the expansion of the bias of the mle) and a simulate method as functions of the model parameters",
##     "gradient of the log-likelihood at the mle", "mle of the dispersion parameter",
##     "expected information matrix evaluated at the mle",
##     "first term in the expansion of the bias of the mle at the mle"),
##     component = list("auxiliary_functions", "score_mle", "dispersion_mle",
##         "expected_information_mle",
##         "bias_mle"), path = "~/Downloads", attempt_rename = FALSE)

