context("test implementation of the simulate auxilliary function")


nsimu <- 100000
tol <- 1e-06

## Binomial
data(endometrial)

eML <- glm(HG ~ NV + PI + EH, family = binomial("probit"), data = endometrial)
theta_mle <- coef(eML)

enriched_eML <- enrich(eML, with = "auxiliary functions")

simu1 <- enriched_eML$auxiliary_functions$simulate(coef(eML), nsim = nsimu, seed = 123)
simu2 <- simulate(eML, nsim = nsimu, seed = 123)

test_that("simulate method and the simulate auxiliary function return the same result with the same seed [glm]", {
    expect_equal(rowMeans(simu1), rowMeans(simu2), tolerance = tol)
})

simu1 <- enriched_eML$auxiliary_functions$simulate(coefficients = c(0.5, 0, 0, 0), nsim = nsimu, seed = 123)
test_that("the simulate auxiliary function does the right thing [glm]", {
    expect_lt(max(abs(rowMeans(simu1) - pnorm(0.5))), 0.01)
})


## Gamma

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))

cML <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)

enriched_cML <- enrich(cML, with = c("mle of dispersion"))
enriched_cML <- enrich(enriched_cML, with = "auxiliary functions")

## Simulation at the ML fit and the ML estimate of dispersion
simu1 <- enriched_cML$auxiliary_functions$simulate(coef(cML),
                                                   dispersion = enriched_cML$dispersion_mle,
                                                   nsim = nsimu, seed = 123)
## The simulate method uses the ML estimate of dispersion only
simu2 <- simulate(cML, nsim = nsimu, seed = 123)

test_that("simulate method and the simulate auxiliary function return the same result with the same seed [glm]", {
    expect_equal(rowMeans(simu1), rowMeans(simu2), tolerance = tol)
})

test_that("simulate returns an error for coefficient vectors with wrong length [glm]", {
    expect_error(enriched_cML$auxiliary_functions$simulate(c(0, 0)))
})

