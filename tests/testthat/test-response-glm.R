context("supplied responses in glm auxiliary functions")

test_that("glm auxiliary functions use a supplied vector response", {
    fit <- glm(counts ~ outcome + treatment, family = poisson(),
               data = data.frame(
                   counts = c(18, 17, 15, 20, 10, 20, 25, 13, 12),
                   outcome = gl(3, 1, 9),
                   treatment = gl(3, 3)))
    aux <- get_auxiliary_functions(fit)
    response <- rev(model.response(model.frame(fit)))
    x <- model.matrix(fit)
    mu <- fitted(fit)

    expect_equal(aux$dmodel(), aux$dmodel(response = model.response(model.frame(fit))))
    expect_equal(aux$pmodel(), aux$pmodel(response = model.response(model.frame(fit))))
    expect_equal(aux$score(), aux$score(response = model.response(model.frame(fit))))
    expect_equal(aux$information(type = "observed"),
                 aux$information(type = "observed",
                                 response = model.response(model.frame(fit))))

    expect_equal(aux$dmodel(response = response),
                 dpois(response, mu), check.attributes = FALSE)
    expect_equal(aux$pmodel(response = response),
                 ppois(response, mu), check.attributes = FALSE)
    expect_equal(aux$score(response = response),
                 drop(crossprod(x, response - mu)),
                 check.attributes = FALSE)
    expect_false(isTRUE(all.equal(aux$score(), aux$score(response = response))))
    expect_equal(aux$information(response = response), aux$information())

    expect_error(aux$score(response = response[-1]), "same length")
    expect_error(aux$dmodel(response = matrix(response, ncol = 1)),
                 "must be a vector")
})

test_that("observed information uses a supplied response", {
    dat <- data.frame(y = c(18, 12, 9, 7, 6, 5, 4, 3),
                      x = seq(-1, 1, length.out = 8))
    fit <- glm(y ~ x, family = Gamma("log"), data = dat)
    aux <- get_auxiliary_functions(fit)
    response <- rev(dat$y)
    dispersion <- summary(fit)$dispersion
    parameters <- c(coef(fit), dispersion)

    supplied_information <- aux$information(
        dispersion = dispersion, type = "observed", response = response)
    numerical_information <- -numDeriv::jacobian(function(parameters) {
        aux$score(coefficients = parameters[-length(parameters)],
                  dispersion = parameters[length(parameters)],
                  response = response)
    }, parameters)

    expect_false(isTRUE(all.equal(
        supplied_information,
        aux$information(dispersion = dispersion, type = "observed"))))
    expect_equal(supplied_information, numerical_information,
                 tolerance = 1e-05, check.attributes = FALSE)
    expect_equal(aux$information(dispersion = dispersion,
                                 response = response),
                 aux$information(dispersion = dispersion))
})

test_that("grouped-binomial responses are initialized once", {
    dat <- data.frame(success = c(2, 4, 3, 7, 5, 8),
                      failure = c(8, 6, 7, 3, 5, 2),
                      x = seq(-1, 1, length.out = 6))
    fit <- glm(cbind(success, failure) ~ x, family = binomial(), data = dat)
    aux <- get_auxiliary_functions(fit)
    response <- cbind(success = c(1, 3, 4, 6, 6, 9),
                      failure = c(9, 7, 6, 4, 4, 1))
    x <- model.matrix(fit)
    mu <- fitted(fit)
    totals <- rowSums(response)

    expect_equal(aux$dmodel(response = response),
                 dbinom(response[, 1], totals, mu),
                 check.attributes = FALSE)
    expect_equal(aux$pmodel(response = response),
                 pbinom(response[, 1], totals, mu),
                 check.attributes = FALSE)
    expect_equal(aux$score(response = response),
                 drop(crossprod(x, response[, 1] - totals * mu)),
                 check.attributes = FALSE)
    expect_equal(aux$information(response = response),
                 crossprod(x * sqrt(totals * mu * (1 - mu))),
                 check.attributes = FALSE)
    expect_equal(aux$qmodel(rep(0.5, nrow(dat))),
                 qbinom(rep(0.5, nrow(dat)), rowSums(dat[1:2]), mu),
                 check.attributes = FALSE)

    expect_error(aux$score(response = response[, 1]), "must be a matrix")
    expect_error(aux$dmodel(response = response[-1, ]), "same dimensions")
    expect_error(aux$qmodel(0.5), "each fitted observation")
})

test_that("factor responses retain the fitted levels", {
    dat <- data.frame(y = factor(rep(c("no", "yes"), 5)),
                      x = seq_len(10))
    fit <- glm(y ~ x, family = binomial(), data = dat)
    aux <- get_auxiliary_functions(fit)
    response <- factor(rep("yes", 10))

    expect_silent(score <- aux$score(response = response))
    expect_length(score, length(coef(fit)))
    expect_equal(aux$dmodel(response = response),
                 dbinom(rep(1, 10), 1, fitted(fit)),
                 check.attributes = FALSE)
    expect_error(aux$score(response = rep("yes", 10)), "must be a factor")
    expect_error(aux$score(response = factor(rep("other", 10))),
                 "levels not present")
})

test_that("responses transformed in the formula are supplied on the model-response scale", {
    dat <- data.frame(y = c(18, 12, 9, 7, 6, 5, 4, 3),
                      x = seq(-1, 1, length.out = 8))
    fit <- glm(log(y) ~ x, family = Gamma("log"), data = dat)
    aux <- get_auxiliary_functions(fit)
    response <- log(dat$y)

    expect_equal(model.response(model.frame(fit)), response,
                 check.attributes = FALSE)
    expect_equal(aux$dmodel(), aux$dmodel(response = response),
                 check.attributes = FALSE)
    expect_equal(aux$pmodel(), aux$pmodel(response = response),
                 check.attributes = FALSE)
    expect_equal(aux$score(), aux$score(response = response))
    expect_equal(aux$information(type = "observed"),
                 aux$information(type = "observed", response = response))
})
