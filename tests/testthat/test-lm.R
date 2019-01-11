context("enrichment of lms")

library("numDeriv")

## see ?longley
longley.x <- data.matrix(longley[, 1:6])
longley.y <- longley[, "Employed"]
## pairs(longley, main = "longley data")
fm1 <- lm(Employed ~ ., data = longley)
fm2 <- glm(Employed ~ ., data = longley)


tol <- 1e-05

test_that("implementation of the scores corresponds to that of the observed information [lm]", {
    enriched_fm1 <- enrich(fm1, with = "auxiliary functions")
    ## MLE of dispersion
    disp <- summary(fm1)$sigma^2 * 9 / 16
    info_appr <- -jacobian(function(coefs) {
        p <- length(coefs)
        enriched_fm1$auxiliary_functions$score(coefs[1:(p-1)],
                                               dispersion = coefs[p])
    }, c(coef(enriched_fm1), disp))
    info_exac <- get_information_function(enriched_fm1)
    expect_equal(info_exac(type = "observed"), info_appr, check.attributes = FALSE)
})


enriched_fm1 <- enrich(fm1, with = "auxiliary functions")
enriched_fm2 <- enrich(fm2, with = "auxiliary functions")
enriched_ufm1 <- enrich(update(fm1, weights = 1:16), with = "auxiliary functions")
enriched_ufm2 <- enrich(update(fm2, weights = 1:16), with = "auxiliary functions")

test_that("impementation of the bias functions corresponds to what glm returns [lm]", {
    expect_equal(enriched_fm1$auxiliary_functions$bias(), enriched_fm2$auxiliary_functions$bias(), tolerance = tol)
    expect_equal(enriched_ufm1$auxiliary_functions$bias(), enriched_ufm2$auxiliary_functions$bias(), tolerance = tol)
})

test_that("simulate returns an error for coefficient vectors with wrong length [lm]", {
    expect_error(enriched_fm1$auxiliary_functions$simulate(c(0, 0, 0)))
})


## ## see ?longley
## longley.x <- data.matrix(longley[, 1:6])
## longley.y <- longley[, "Employed"]
## fm1 <- glm(Employed ~ ., data = longley)
## fm2 <- lm(Employed ~ ., data = longley)
## ufm1 <- update(fm1, weights = 1:16)
## ufm2 <- update(fm2, weights = 1:16)

## ## The variance estimate for the unweighted fit is
## summary(fm1)$dispersion
## summary(fm2)$sigma^2
## ## while for the weighted one is
## summary(ufm1)$dispersion
## summary(ufm2)$sigma^2

## ## The variance estimate for the weighted fit is simply
## sum(residuals(ufm1)^2)/ufm1$df.residual
## ## or
## sum((ufm1$y - fitted(ufm1))^2 * ufm1$prior.weights)/ufm1$df.residual
## ## because
## all.equal(residuals(ufm1), (ufm1$y - fitted(ufm1)) * sqrt(ufm1$prior.weights))



## enriched_fm1 <- enrich(fm1, with = "auxiliary functions")
## disp <- summary(fm1)$dispersion

## info_appr <- -jacobian(function(coefs) {
##     p <- length(coefs)
##     enriched_fm1$auxiliary_functions$score(coefs[1:(p-1)],
##                                            dispersion = coefs[p])
## }, c(coef(enriched_fm1), disp))
## info_exac <- get_information_function(enriched_fm1)
## solve(info_exac(type = "observed")) -  solve(info_appr)
