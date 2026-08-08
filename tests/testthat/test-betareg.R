context("enrichment of betareg objects")


library("betareg")
library("numDeriv")

## Section 4 from Ferrari and Cribari-Neto (2004)
data("GasolineYield", package = "betareg")
data("FoodExpenditure", package = "betareg")
tol <- 1e-08

test_that("extended-support beta regressions are rejected", {
    fit <- betareg(yield ~ batch + temp, data = GasolineYield)

    for (dist in c("xbeta", "xbetax")) {
        fit$dist <- dist
        expect_error(get_auxiliary_functions(fit), 'dist = "beta"',
                     fixed = TRUE)
    }
})

for (link in c("logit", "probit", "cloglog", "cauchit")) {
    gy <- betareg(yield ~ batch + temp | temp, data = GasolineYield, link = link)

    gy_enriched <- enrich(gy)

    test_that("implementation of the expected information is correct [betareg]",
    {
        expect_equal(solve(get_information_function(gy)()), vcov(gy), tolerance = tol, check.attributes = FALSE)
    })


    test_that("implementation of the scores corresponds to that of the observed information [betareg]",
    {
        info_appr <- -jacobian(function(coefs) {
                          gy_enriched$auxiliary_functions$score(coefs)
                      },
                      x = c(coef(gy)))
        info_exac <- gy_enriched$auxiliary_functions$information(coef(gy_enriched),
                                                                 type = "observed")
        expect_equal(solve(info_appr), solve(info_exac), tolerance = 1e-03, check.attributes = FALSE)
    })

    test_that("score and information are derivatives of the dmodel log-likelihood [betareg]", {
        auxiliary <- gy_enriched$auxiliary_functions
        coefficients <- c(coef(gy))
        coefficients[1] <- coefficients[1] + 0.1
        case_weights <- gy$weights
        if (is.null(case_weights)) {
            case_weights <- rep.int(1, nobs(gy))
        }
        log_likelihood <- function(coefficients) {
            sum(case_weights * auxiliary$dmodel(
                coefficients = coefficients, log = TRUE))
        }

        expect_equal(auxiliary$score(coefficients),
                     grad(log_likelihood, coefficients),
                     tolerance = 1e-05, check.attributes = FALSE)
        expect_equal(auxiliary$information(coefficients, type = "observed"),
                     -hessian(log_likelihood, coefficients),
                     tolerance = 1e-05, check.attributes = FALSE)
    })

    test_that("implementation of the bias function corresponds to what is computed internally [betareg]", {
        b_betareg <- update(gy, type = "BC")$bias
        b_enrichwith <- get_bias_function(gy)(coef(gy))
        expect_equal(b_betareg, b_enrichwith)
    })

    test_that("score contributions from enrichwith are the same to those computed internally [betareg]", {
        scores <- gy_enriched$auxiliary_functions$score
        expect_equal(scores(contributions = TRUE), gy$grad)
    })

    test_that("score contributions add up to scores [betareg]",
    {
        scores <- gy_enriched$auxiliary_functions$score
        expect_equal(colSums(scores(contributions = TRUE)), scores())
    })

    test_that("simulation function works [betareg]",
    {
        bsimulate <- gy_enriched$auxiliary_functions$simulate
        expect_error(bsimulate(c(0, 1, 2)))
        expect_true(isTRUE(identical(nrow(bsimulate(nsim = 100)), nrow(gy$model))))
        expect_true(isTRUE(identical(ncol(bsimulate(nsim = 99)), as.integer(99))))
    })
}
