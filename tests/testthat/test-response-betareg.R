context("supplied responses in betareg auxiliary functions")

test_that("betareg score and information use a supplied response", {
    data("GasolineYield", package = "betareg")
    fit <- betareg::betareg(yield ~ batch + temp | temp,
                           data = GasolineYield)
    aux <- get_auxiliary_functions(fit)
    fitted_response <- fit$y
    response <- rev(fitted_response)
    coefficients <- coef(fit, model = "full")

    expect_equal(aux$score(), aux$score(response = fitted_response))
    expect_equal(aux$information(type = "observed"),
                 aux$information(type = "observed",
                                 response = fitted_response))
    expect_false(isTRUE(all.equal(aux$score(),
                                  aux$score(response = response))))
    expect_equal(aux$information(response = response), aux$information())
    expect_false(isTRUE(all.equal(
        aux$information(type = "observed", response = response),
        aux$information(type = "observed"))))

    numerical_information <- -numDeriv::jacobian(function(coefficients) {
        aux$score(coefficients = coefficients, response = response)
    }, coefficients)
    expect_equal(aux$information(coefficients = coefficients,
                                 type = "observed", response = response),
                 numerical_information, tolerance = 1e-05,
                 check.attributes = FALSE)

    expect_error(aux$score(response = response[-1]), "same length")
    expect_error(aux$score(response = matrix(response, ncol = 1)),
                 "numeric vector")
    expect_error(aux$information(response = replace(response, 1, 0)),
                 "in \\(0, 1\\)")
    expect_error(aux$information(response = replace(response, 1, NA_real_)),
                 "finite and in")
})

test_that("betareg dmodel, pmodel, and qmodel use the fitted design", {
    data("GasolineYield", package = "betareg")
    fit <- betareg::betareg(yield ~ batch + temp | temp,
                           data = GasolineYield)
    aux <- get_auxiliary_functions(fit)
    response <- rev(fit$y)
    coefficients <- coef(fit, model = "full")
    coefficients[1] <- coefficients[1] + 0.1
    k <- ncol(model.matrix(fit, model = "mean"))
    beta <- coefficients[seq_len(k)]
    gamma <- coefficients[k + seq_len(
        ncol(model.matrix(fit, model = "precision")))]
    eta <- drop(model.matrix(fit, model = "mean") %*% beta)
    phi_eta <- drop(model.matrix(fit, model = "precision") %*% gamma)
    mu <- fit$link$mean$linkinv(eta)
    phi <- fit$link$precision$linkinv(phi_eta)
    shape1 <- mu * phi
    shape2 <- (1 - mu) * phi
    probabilities <- seq(0.05, 0.95, length.out = nobs(fit))

    expect_equal(aux$dmodel(), aux$dmodel(response = fit$y))
    expect_equal(aux$pmodel(), aux$pmodel(response = fit$y))
    expect_equal(aux$dmodel(response, coefficients = coefficients),
                 dbeta(response, shape1, shape2),
                 check.attributes = FALSE)
    expect_equal(aux$pmodel(response, coefficients = coefficients),
                 pbeta(response, shape1, shape2),
                 check.attributes = FALSE)
    expect_equal(aux$pmodel(response, coefficients = coefficients,
                            lower.tail = FALSE, log.p = TRUE),
                 pbeta(response, shape1, shape2,
                       lower.tail = FALSE, log.p = TRUE),
                 check.attributes = FALSE)
    expect_equal(aux$qmodel(probabilities, coefficients = coefficients),
                 qbeta(probabilities, shape1, shape2),
                 check.attributes = FALSE)
    expect_equal(aux$qmodel(log(probabilities), coefficients = coefficients,
                            log.p = TRUE),
                 qbeta(log(probabilities), shape1, shape2, log.p = TRUE),
                 check.attributes = FALSE)
    expect_equal(get_dmodel_function(fit)(response), aux$dmodel(response))
    expect_equal(get_pmodel_function(fit)(response), aux$pmodel(response))
    expect_equal(get_qmodel_function(fit)(probabilities),
                 aux$qmodel(probabilities))
    expect_error(aux$qmodel(probabilities[-1]), "each fitted observation")
})

test_that("transformed betareg responses are supplied on the model-response scale", {
    data("GasolineYield", package = "betareg")
    fit <- betareg::betareg(I(yield / 2) ~ batch + temp,
                           data = GasolineYield)
    aux <- get_auxiliary_functions(fit)
    response <- GasolineYield$yield / 2

    expect_equal(aux$score(), aux$score(response = response))
    expect_equal(aux$information(type = "observed"),
                 aux$information(type = "observed", response = response))
})
