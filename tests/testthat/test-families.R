context("implementation of derivatives of variance function with respect to mu (comparison with numerical derivatives)")

library("numDeriv")

families <- list(binomial(),
                 poisson(),
                 Gamma(),
                 inverse.gaussian(),
                 gaussian(),
                 quasi(variance = "constant"),
                 quasi(variance = "mu(1-mu)"),
                 quasi(variance = "mu"),
                 quasi(variance = "mu^2"),
                 quasi(variance = "mu^3"),
                 quasibinomial(),
                 quasipoisson())

tol <- 0.00001 # sqrt(.Machine$double.eps)
epsilon <- 1e-10
len <- 100


## Test variance and derivatives
for (fam in families) {

    ## A range for mu
    mus <- seq(-10, 10, length = len)
    cfam <- enrich(fam)

    test_that(paste("d1variance is correctly implemented for", fam$family), {
        expect_equal(grad(cfam$variance, mus), cfam$d1variance(mus), tolerance = tol)
    })

    test_that(paste("d2variance is correctly implemented for", fam$family), {
        expect_equal(grad(cfam$d1variance, mus), cfam$d2variance(mus), tolerance = tol)
    })
}


## Test afun and derivatives
for (fam in families) {
    if (fam$family %in% c("poisson", "binomial", "quasi", "quasibinomial", "quasipoisson")) {
        next
    }
    ## A range for zeta
    zetas <- seq(-3, -0.01, length = len)
    cfam <- enrich(fam)

    test_that(paste("d2afun is correctly implemented for", fam$family), {
        expect_equal(grad(cfam$d1afun, zetas), cfam$d2afun(zetas), tolerance = tol)
    })

    test_that(paste("d3afun is correctly implemented for", fam$family), {
        expect_equal(grad(cfam$d2afun, zetas), cfam$d3afun(zetas), tolerance = tol)
    })

}
