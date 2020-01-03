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

## Create a test data set
test_data <- model1$data[1:8, ]
test_data$grahami <- c(1, 1, 1, 1, 0, 0, 0, 0)
test_data$opalinus <- c(0, 0, 0, 0, 1, 1, 1, 1)
test_data$species <- factor(rep(c("grahami", "opalinus"), each = 4), levels = c("opalinus", "grahami"))

test_that("dmodel returns the same results for various equivalent representations of the data for logistic regression", {
    expect_equal(enrich(model1)$auxiliary_functions$dmodel(test_data),
                 enrich(model3)$auxiliary_functions$dmodel(test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
    expect_equal(enrich(model1)$auxiliary_functions$dmodel(test_data),
                 enrich(model2)$auxiliary_functions$dmodel(test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
    expect_equal(enrich(model2)$auxiliary_functions$dmodel(test_data),
                 enrich(model3)$auxiliary_functions$dmodel(test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
})


test_that("pmodel returns the same results for various equivalent representations of the data for logistic regression", {
    expect_equal(enrich(model1)$auxiliary_functions$pmodel(test_data),
                 enrich(model3)$auxiliary_functions$pmodel(test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
    expect_equal(enrich(model1)$auxiliary_functions$pmodel(test_data),
                 enrich(model2)$auxiliary_functions$pmodel(test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
    expect_equal(enrich(model2)$auxiliary_functions$pmodel(test_data),
                 enrich(model3)$auxiliary_functions$pmodel(test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
})


ps <- seq(0, 1, length = nrow(test_data))
test_that("qmodel returns the same results for various equivalent representations of the data for logistic regression", {
    expect_equal(enrich(model1)$auxiliary_functions$qmodel(ps, test_data),
                 enrich(model3)$auxiliary_functions$qmodel(ps, test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
    expect_equal(enrich(model1)$auxiliary_functions$qmodel(ps, test_data),
                 enrich(model2)$auxiliary_functions$qmodel(ps, test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
    expect_equal(enrich(model2)$auxiliary_functions$qmodel(ps, test_data),
                 enrich(model3)$auxiliary_functions$qmodel(ps, test_data),
                 check.attributes = FALSE,
                 tol = tolerance)
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
