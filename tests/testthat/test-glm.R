context("maximum likelihood estimate of gamma dispersion")

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)

tol <- 1e-06
test_that("ML estimate of gamma dispersion from enrichwith is numerically the same to that from MASS::gamma.dispersion", {
    enriched_mod <- enrich(mod, with = "MLE of dispersion")
    expect_equal(enriched_mod$dispersion.mle,
                 MASS::gamma.dispersion(mod),
                 tolerance = tol)
    })
