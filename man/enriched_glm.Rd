% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enriched_glm.R
\name{enriched_glm}
\alias{enriched_glm}
\title{Fitting generalized linear models enriched with useful components}
\usage{
enriched_glm(formula, family = gaussian, ...)
}
\arguments{
\item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

\item{family}{a description of the error distribution and link
    function to be used in the model.  For \code{glm} this can be a
    character string naming a family function, a family function or the
    result of a call to a family function.  For \code{glm.fit} only the
    third option is supported.  (See \code{\link{family}} for details of
    family functions.)}

\item{...}{other arguments passed to \code{\link{glm}}}
}
\value{
An object of class \code{enriched_glm} that contains all the
components of a \code{\link{glm}} object, along with a set of
auxiliary functions (score function, information matrix, a simulate
method, first term in the expansion of the bias of the maximum
likelihood estimator, and dmodel, pmodel, qmodel), the maximum
likelihood estimate of the dispersion parameter, the expected or
observed information at the maximum likelihood estimator, and the
first term in the expansion of the bias of the maximum likelihood
estimator.

See \code{\link{enrich.glm}} for more details and links for the
auxiliary functions.
}
\description{
\code{\link{enriched_glm}} fits generalized linear models using
\code{\link{glm}} and then enriches the resulting object with all
enrichment options.
}
\details{
\code{enriched_glm} has the same interface as \code{\link{glm}}
}
\examples{
\dontrun{
# A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
   u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
   time = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
   lot = factor(c(rep(1, 9), rep(2, 9))))

# Fit a generalized linear model
cML <- enriched_glm(time ~ lot*log(u), data = clotting, family = Gamma("log"))

# Evaluate the densities at the data points in clotting at the
# maximum likelihood estimates
cML_dmodel <- get_dmodel_function(cML) # same as cML$auxiliary_functions$dmodel
cML_dmodel()

# Evaluate the densities at supplied data points
new_data <- data.frame(u = c(15:17, 15:17),
                       time = c(30:32, 15:17),
                       lot = factor(c(1, 1, 1, 2, 2, 2)))
cML_dmodel(data = new_data)

# Get pmodel and qmodel function
cML_pmodel <- get_pmodel_function(cML) # same as cML$auxiliary_functions$pmodel
cML_qmodel <- get_qmodel_function(cML) # same as cML$auxiliary_functions$qmodel

# The following should return c(30:32, 15:17)
probs <- cML_pmodel(data = new_data)
cML_qmodel(probs, data = new_data)

# Evaluate the observed information matrix at the MLE
cML_info <- get_information_function(cML)
cML_info(type = "observed")

# Wald tests based on the observed information at the
# moment based esimator of the dispersion
dispersion <- summary(cML)$dispersion
cML_vcov_observed <- solve(cML_info(dispersion = dispersion, type = "observed"))
lmtest::coeftest(cML, vcov = cML_vcov_observed)

# Wald tests based on the expected information at the
# moment based esimator of the dispersion
cML_vcov_expected <- solve(cML_info(dispersion = dispersion, type = "expected"))
lmtest::coeftest(cML, vcov = cML_vcov_expected)
# Same statistics as coef(summary(cML))[, "t value"]

}


}
