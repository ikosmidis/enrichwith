context("enrichment of glms")

library("numDeriv")
library("MASS")

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod1 <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)

tol <- 1e-05
test_that("implementation of the scores corresponds to that of the observed information",
{
    enriched_mod1 <- enrich(mod1, with = "auxiliary functions")
    disp <- summary(enriched_mod1)$disp
    info_appr <- -jacobian(function(coefs) {
        p <- length(coefs)
        enriched_mod1$auxiliary_functions$score(coefs[1:(p-1)],
                                          dispersion = coefs[p])
    }, c(coef(enriched_mod1), disp))
    info_exac <- enriched_mod1$auxiliary_functions$information(coef(enriched_mod1),
                                                         dispersion = disp,
                                                         type = "observed")
    expect_equal(sum(abs(info_appr - info_exac)), 0, tolerance = 0.01)
})

test_that("ML estimate of gamma dispersion from enrichwith is numerically the same to that from MASS::gamma.dispersion", {
    enriched_mod1 <- enrich(mod1, with = "mle of dispersion")
    expect_equal(unname(enriched_mod1$dispersion_mle),
                 gamma.dispersion(mod1),
                 tolerance = tol)
})

test_that("standard errors from expected information matrix match those from summary", {
    dispersion_pearson <- summary(mod1)$dispersion
    enriched_mod1 <- enrich(mod1, with = "expected information", dispersion = dispersion_pearson)
    stderrs1 <- sqrt(diag(solve(enriched_mod1$expected_information)))
    stderrs2 <- coef(summary(mod1))[, 2]
    expect_equal(stderrs1[-length(stderrs1)], stderrs2)
})


## A Poisson example
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
mod2 <- glm(counts ~ outcome + treatment, family = poisson())

test_that("expected information matrix from enrich with is equal to that coming from the qr decomposition component of the object", {
    enriched_mod2 <- enrich(mod2, with = "auxiliary functions")
    expect_equal(enriched_mod2$auxiliary_functions$information(coef(mod2), type = "expected"),
                 crossprod(qr.R(mod2$qr)),
                 tolerance = 200)
})



