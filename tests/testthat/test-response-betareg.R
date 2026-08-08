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
