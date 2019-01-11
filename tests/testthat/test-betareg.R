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
                      gy_enriched$auxiliary_functions$score(coefs)
                  },
                  x = c(coef(gy)),
                  method.args = list(eps = 1e-7, d = 0.01,
                                     zero.tol=sqrt(.Machine$double.eps/7e-7), r = 10, v = 5))
    info_exac <- gy_enriched$auxiliary_functions$information(coef(gy_enriched),
                                                               type = "observed")
    expect_equal(solve(info_appr), solve(info_exac), tolerance = 1e-03, check.attributes = FALSE)
})

test_that("score contributions add up to scores",
{
    scores <- gy_enriched$auxiliary_functions$score
    expect_equal(colSums(scores(contributions = TRUE)), scores())
})
