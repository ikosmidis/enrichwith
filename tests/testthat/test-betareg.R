context("enrichment of betareg objects")


library("betareg")


## Section 4 from Ferrari and Cribari-Neto (2004)
data("GasolineYield", package = "betareg")
data("FoodExpenditure", package = "betareg")

## Table 1
gy <- betareg(yield ~ batch + temp | temp, data = GasolineYield)

gy_enriched <- enrich(gy)

tol <- 1e-08
test_that("implementation of the expected information is correct",
{
    expect_equal(solve(get_information_function(gy)()), vcov(gy), tolerance = tol, check.attributes = FALSE)
})

test_that("implementation of the scores corresponds to that of the observed information",
{
    info_appr <- -jacobian(function(coefs) {
        p <- length(coefs)
        gy_enriched$auxiliary_functions$score(coefs)
    }, c(coef(gy_enriched)))
    info_exac <- gy_enriched$auxiliary_functions$information(coef(gy_enriched),
                                                               type = "observed")
    expect_equal(info_appr, info_exac, tolerance = 0.01)
})

