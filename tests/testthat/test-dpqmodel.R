context("implementation of model auxiliary functions")

## Tolerance for comparisons
tolerance <- 1e-05

data("lizards", package = "brglm")

model1 <- glm(formula = grahami/(grahami + opalinus) ~ height + diameter + light + time, family = binomial(logit), weights = (grahami + opalinus), data = lizards)
model2 <- glm(formula = cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(logit), data = lizards)

lizards_grahami <- lizards[, c("grahami", "height", "diameter", "light", "time")]
lizards_grahami <- lizards_grahami[rep(seq.int(nrow(lizards_grahami)), lizards_grahami$grahami), ]
lizards_grahami$species <- "grahami"
lizards_grahami$grahami <- NULL
lizards_opalinus <- lizards[, c("opalinus", "height", "diameter", "light", "time")]
lizards_opalinus <- lizards_opalinus[rep(seq.int(nrow(lizards_opalinus)), lizards_opalinus$opalinus), ]
lizards_opalinus$species <- "opalinus"
lizards_opalinus$opalinus <- NULL
lizards1 <- rbind(lizards_grahami, lizards_opalinus)
lizards1$species <- factor(lizards1$species, levels = c("opalinus", "grahami"))

model3 <- glm(formula = species ~ height + diameter + light + time, family = binomial(logit), data = lizards1)

test_that("simulate and get_simulate_function return the same variates for various equivalent representations of the data for logistic regression", {
    expect_identical(simulate(model1, seed = 123)[, 1],
                     get_simulate_function(model1)(seed = 123)[, 1])
    expect_identical(simulate(model2, seed = 123)[, 1],
                     get_simulate_function(model2)(seed = 123)[, 1])
    expect_identical(simulate(model3, seed = 123)[, 1],
                     get_simulate_function(model3)(seed = 123)[, 1])
})

probs <- 1:10 / 11
test_that("qmodel returns the same output across data representations", {
    for (pr in probs) {
        expect_equal(aux1$qmodel(rep(pr, nrow(lizards))), aux2$qmodel(rep(pr, nrow(lizards))))
    }
})

aux1 <- get_auxiliary_functions(model1)
aux2 <- get_auxiliary_functions(model2)
aux3 <- get_auxiliary_functions(model3)
tots <- lizards$grahami + lizards$opalinus
d1 <- simulate(model1, seed = 123)[, 1]
d2 <- simulate(model2, seed = 123)[, 1]

test_that("d/pmodel returns the same results for various equivalent representations of the data for logistic regression", {
    expect_equal(aux1$dmodel(d2[, 1] / rowSums(d2)) , aux1$dmodel(d1), check.attributes = FALSE)
    expect_equal(aux2$dmodel(cbind(d1, 1- d1) * tots), aux2$dmodel(d2), check.attributes = FALSE)
    expect_equal(aux1$pmodel(d2[, 1] / rowSums(d2)) , aux1$pmodel(d1), check.attributes = FALSE)
    expect_equal(aux2$pmodel(cbind(d1, 1- d1) * tots), aux2$pmodel(d2), check.attributes = FALSE)
})


library("numDeriv")
library("MASS")

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod1 <- glm(conc ~ lot*log(u), data = clotting, family = inverse.gaussian)

test_that("d/p/r/qmodel works for inverse gaussian regression", {
    simulate_ig <- get_simulate_function(mod1)
    set.seed(123)
    s1 <- simulate_ig(coefficients = coef(mod1), dispersion = summary(mod1)$dispersion)
    set.seed(123)
    s2 <- simulate(mod1)
    expect_equal(s1, s2, check.attributes = FALSE)

    disp <- enrich(mod1)$dispersion

    d1 <- SuppDists::dinvGauss(clotting$conc, fitted.values(mod1), lambda = 1/disp)
    p1 <- SuppDists::pinvGauss(clotting$conc, fitted.values(mod1), lambda = 1/disp)
    q1 <- SuppDists::qinvGauss(rep(0.2, 18), fitted.values(mod1), lambda = 1/disp)

    d2 <- get_dmodel_function(mod1)()
    p2 <- get_pmodel_function(mod1)()
    q2 <- get_qmodel_function(mod1)(rep(0.2, 18))

    expect_equal(d1, d2, check.attributes = FALSE)
    expect_equal(q1, q2, check.attributes = FALSE)
    expect_equal(p1, p2, check.attributes = FALSE)


})
