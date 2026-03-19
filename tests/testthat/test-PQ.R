context("Implementation of P/Q matrices")

library("numDeriv")
library("MASS")

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))

test_that("bias implementation matches manual implementation through P and Q [Gamma]", {
    mod1 <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)
    mod1e <- enrich(mod1, with = "auxiliary functions")
    v <- solve(mod1e$auxiliary_functions$information())
    P <- mod1e$auxiliary_functions$Pmat()
    Q <- mod1e$auxiliary_functions$Qmat()
    b0 <- mod1e$auxiliary_functions$bias()
    expect_equal(b0, - drop(v %*% sapply(seq.int(length(coef(mod1)) + 1), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)

    coefs <- coef(mod1e)
    coefs <- coefs * 10
    disp <- 0.5
    v <- solve(mod1e$auxiliary_functions$information(coefs, disp))
    P <- mod1e$auxiliary_functions$Pmat(coefs, disp)
    Q <- mod1e$auxiliary_functions$Qmat(coefs, disp)
    b0 <- mod1e$auxiliary_functions$bias(coefs, disp)
    expect_equal(b0, - drop(v %*% sapply(seq.int(length(coef(mod1)) + 1), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)

})


test_that("enriching a brglmFit object returns same reults as enriching a glm object and evaluating at RB estimates", {
    mML <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)
    mRB <- glm(conc ~ lot*log(u), data = clotting, family = Gamma, method = brglm2::brglmFit)
    eML <- enrich(mML, with = "auxiliary functions")
    eRB <- enrich(mRB, with = "auxiliary functions")
    expect_equal(eRB$auxiliary_functions$information(dispersion = coef(mRB, model = "dispersion")),
                 eML$auxiliary_functions$information(coef(mRB), coef(mRB, model = "dispersion")))
    expect_equal(eRB$auxiliary_functions$Pmat(dispersion = coef(mRB, model = "dispersion")),
                 eML$auxiliary_functions$Pmat(coef(mRB), coef(mRB, model = "dispersion")))
    expect_equal(eRB$auxiliary_functions$Qmat(dispersion = coef(mRB, model = "dispersion")),
                 eML$auxiliary_functions$Qmat(coef(mRB), coef(mRB, model = "dispersion")))
})


## A binomial examples
data("lizards", package = "brglm2")

test_that("bias implementation matches manual implementation through P and Q [binomial(cauchit)]", {
    lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
                     light + time, family = binomial(cauchit), data = lizards)
    mod1e <- enrich(lizardsML, with = "auxiliary functions")
    v <- solve(mod1e$auxiliary_functions$information())
    P <- mod1e$auxiliary_functions$Pmat()
    Q <- mod1e$auxiliary_functions$Qmat()
    b0 <- mod1e$auxiliary_functions$bias()
    expect_equal(b0, - drop(v %*% sapply(seq.int(length(coef(mod1e))), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)
    coefs <- coef(mod1e)
    coefs <- coefs * 10
    disp <- 0.5
    v <- solve(mod1e$auxiliary_functions$information(coefs, disp))
    P <- mod1e$auxiliary_functions$Pmat(coefs, disp)
    Q <- mod1e$auxiliary_functions$Qmat(coefs, disp)
    b0 <- mod1e$auxiliary_functions$bias(coefs, disp)
    expect_equal(b0, - drop(v %*% sapply(seq.int(length(coef(mod1e))), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)
})

test_that("bias implementation matches manual implementation through P and Q [binomial(logit)]", {
    lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
                         light + time, family = binomial(logit), data = lizards)
    mod1e <- enrich(lizardsML, with = "auxiliary functions")
    coefs <- coef(mod1e)
    nvars <- length(coefs)
    cnams <- names(coefs)
    v <- solve(mod1e$auxiliary_functions$information())
    P <- mod1e$auxiliary_functions$Pmat()
    Q <- mod1e$auxiliary_functions$Qmat()
    b0 <- mod1e$auxiliary_functions$bias()
    expect_equal(b0, - drop(v %*% sapply(seq.int(nvars), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)
    coefs <- coefs * 10
    disp <- 0.5
    v <- solve(mod1e$auxiliary_functions$information(coefs, disp))
    P <- mod1e$auxiliary_functions$Pmat(coefs, disp)
    Q <- mod1e$auxiliary_functions$Qmat(coefs, disp)
    b0 <- mod1e$auxiliary_functions$bias(coefs, disp)
    expect_equal(b0, - drop(v %*% sapply(seq.int(nvars), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)
    for (k in seq.int(nvars))
        expect_equal(Q[[k]], matrix(0, nvars, nvars, dimnames = list(cnams, cnams)))
})



## A Poisson example
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

test_that("bias implementation matches manual implementation through P and Q [poisson(log)]", {
    mod2 <- glm(counts ~ outcome + treatment, family = poisson())
    enriched_mod2 <- enrich(mod2, with = "auxiliary functions")
    mod2e <- enrich(mod2, with = "auxiliary functions")
    coefs <- coef(mod2e)
    nvars <- length(coefs)
    cnams <- names(coefs)
    v <- solve(mod2e$auxiliary_functions$information())
    P <- mod2e$auxiliary_functions$Pmat()
    Q <- mod2e$auxiliary_functions$Qmat()
    b0 <- mod2e$auxiliary_functions$bias()
    expect_equal(b0, - drop(v %*% sapply(seq.int(nvars), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)
    for (k in seq.int(nvars))
        expect_equal(Q[[k]], matrix(0, nvars, nvars, dimnames = list(cnams, cnams)))
})

test_that("enriching a brglmFit object returns same reults as enriching a glm object and evaluating at RB estimates", {
    mML <- glm(counts ~ outcome + treatment, family = poisson())
    mRB <- glm(counts ~ outcome + treatment, family = poisson(), method = brglm2::brglmFit)
    eML <- enrich(mML, with = "auxiliary functions")
    eRB <- enrich(mRB, with = "auxiliary functions")
    expect_equal(eRB$auxiliary_functions$information(),
                 eML$auxiliary_functions$information(coef(mRB)))
    expect_equal(eRB$auxiliary_functions$Pmat(),
                 eML$auxiliary_functions$Pmat(coef(mRB)))
    expect_equal(eRB$auxiliary_functions$Qmat(),
                 eML$auxiliary_functions$Qmat(coef(mRB)))
})

test_that("bias implementation matches manual implementation through P and Q [poisson(sqrt)]", {
    mod2 <- glm(counts ~ outcome + treatment, family = poisson("sqrt"))
    enriched_mod2 <- enrich(mod2, with = "auxiliary functions")
    mod2e <- enrich(mod2, with = "auxiliary functions")
    coefs <- coef(mod2e)
    nvars <- length(coefs)
    cnams <- names(coefs)
    v <- solve(mod2e$auxiliary_functions$information())
    P <- mod2e$auxiliary_functions$Pmat()
    Q <- mod2e$auxiliary_functions$Qmat()
    b0 <- mod2e$auxiliary_functions$bias()
    expect_equal(b0, - drop(v %*% sapply(seq.int(nvars), function(t) sum(v * (P[[t]] + Q[[t]])) / 2)),
                 check.attributes = FALSE)
})


data("coalition", package = "brglm2")

test_that("bias implementation matches manual implementation through P and Q [Gamma - coalition]", {
    mod1 <- glm(duration ~ fract + I(2 * fract) + numst2, family = Gamma, data = coalition)
    mod1e <- enrich(mod1, with = "auxiliary functions")
    i <- mod1e$auxiliary_functions$information()
    P <- mod1e$auxiliary_functions$Pmat()
    Q <- mod1e$auxiliary_functions$Qmat()
    b0 <- mod1e$auxiliary_functions$bias()
    na_coefs <- is.na(coef(mod1e))
    v <- i
    i <- i[!na_coefs, !na_coefs]
    v[!na_coefs, !na_coefs] <- solve(i)
    expect_true(all(sapply(P, function(mat) all(is.na(mat[which(na_coefs), ])))))
    expect_true(all(sapply(P, function(mat) all(is.na(mat[, which(na_coefs)])))))
    expect_true(all(sapply(Q, function(mat) all(is.na(mat[which(na_coefs), ])))))
    expect_true(all(sapply(Q, function(mat) all(is.na(mat[, which(na_coefs)])))))
    b1 <- - colSums(v * sapply(seq.int(length(coef(mod1)) + 1), function(t) sum(v * (P[[t]] + Q[[t]]), na.rm = TRUE) / 2), na.rm = TRUE)
    b1[na_coefs] <- NA
    expect_equal(b0, b1, check.attributes = FALSE)
})
