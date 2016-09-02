Introduction
============

The **enrichwith** package provides the `enrich` method to enrich
list-like R objects with new, relevant components. The resulting objects
preserve their class, so all methods associated with them still apply.

This short vignette demonstrates some of the available enrichment
options for `glm` objects in
[**enrichwith**](https://github.com/ikosmidis/enrichwith). In
particular, the enriched `glm` objects are used here to implement a
quasi Fisher scoring procedure for computing reduced-bias estimates in
[generalized linear
models](https://en.wikipedia.org/wiki/Generalized_linear_model) (see
Kosmidis and Firth 2010 for a parallel quasi
[Newton-Raphson](https://en.wikipedia.org/wiki/Newton%27s_method)
procedure).

------------------------------------------------------------------------

Endometrial cancer data
=======================

Heinze and Schemper (2002) used a logistic regression model to analyse
data from a study on endometrial cancer. Agresti (2015 Section 5.7)
provide details on the data set. Below, we fit a probit regression model
with the same linear predictor as the logistic regression model in
Heinze and Schemper (2002).

    # Get the data from the online supplmementary material of Agresti (2015)
    endometrial_url <- url("http://www.stat.ufl.edu/~aa/glm/data/Endometrial.dat")
    endometrial <- read.table(endometrial_url, header = TRUE)
    modML <- glm(HG ~ NV + PI + EH, family = binomial("probit"), data = endometrial)
    theta_mle <- coef(modML)
    summary(modML)

    ## 
    ## Call:
    ## glm(formula = HG ~ NV + PI + EH, family = binomial("probit"), 
    ##     data = endometrial)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.47007  -0.67917  -0.32978   0.00008   2.74898  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   2.18093    0.85732   2.544 0.010963 *  
    ## NV            5.80468  402.23641   0.014 0.988486    
    ## PI           -0.01886    0.02360  -0.799 0.424066    
    ## EH           -1.52576    0.43308  -3.523 0.000427 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 104.90  on 78  degrees of freedom
    ## Residual deviance:  56.47  on 75  degrees of freedom
    ## AIC: 64.47
    ## 
    ## Number of Fisher Scoring iterations: 17

As is the case for the logistic regression model in Heinze and Schemper
(2002), the maximum likelihood estimate of the parameter for `NV` is
actually infinite. The reported, apparently finite value is merely due
to false convergence of the iterative estimation procedure. The same is
true for the estimated standard error, and, hence the value 0.014 for
the *z*-statistic cannot be trusted for inference on the size of the
effect for `NV`.

For generalized linear models with categorical response models like the
above, the bias reduction method in Firth (1993) has been found to
result in finite estimates even when the maximum likelihood ones are
infinite (see, Heinze and Schemper 2002, for logistic regressions;
Kosmidis and Firth 2011, for multinomial regressions; Kosmidis 2014, for
cumulative link models).

The [**brglm**](https://cran.r-project.org/web/packages/brglm/brglm.pdf)
R package implements one of the variants of the bias reduction method in
Firth (1993) for binomial response generalized linear models and using
iterative maximum likelihood fits on binomial pseudo-data (see, Kosmidis
2007, Chapter 5, for details). The reduced-bias estimates for the probit
regression on the endometrial data can be computed as follows.

    library("brglm")
    modBR <- brglm(HG ~ NV + PI + EH, family = binomial("probit"), data = endometrial)
    theta_brglm <- coef(modBR)
    summary(modBR)

    ## 
    ## Call:
    ## brglm(formula = HG ~ NV + PI + EH, family = binomial("probit"), 
    ##     data = endometrial)
    ## 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  1.91460    0.78877   2.427 0.015210 *  
    ## NV           1.65892    0.74730   2.220 0.026427 *  
    ## PI          -0.01520    0.02089  -0.728 0.466793    
    ## EH          -1.37988    0.40329  -3.422 0.000623 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 93.983  on 78  degrees of freedom
    ## Residual deviance: 57.587  on 75  degrees of freedom
    ## AIC:  65.587

The *z*-statistic for `NV` has value 2.22 when based on the reduced-bias
estimates, providing some evidence for the existence of an effect.

In the following, we use **enrichwith** to implement two variants of the
bias reduction method via a unifying quasi Fisher scoring estimation
procedure.

------------------------------------------------------------------------

Bias reduction in a nutshell
============================

Consider a parametric [statistical
model](https://en.wikipedia.org/wiki/Statistical_model) 𝒫<sub>*θ*</sub>
with unknown parameter *θ* ∈ ℜ<sup>*p*</sup> and the iteration
*θ*<sup>(*k* + 1)</sup> := *θ*<sup>(*k*)</sup> + {*i*(*θ*<sup>(*k*)</sup>)}<sup>−1</sup>*s*(*θ*<sup>(*k*)</sup>)−*c*(*θ*<sup>(*k*)</sup>)*b*(*θ*<sup>(*k*)</sup>)
 where *θ*<sup>(*k*)</sup> is the value of *θ* at the *k*th iteration,
*s*(*θ*) is the gradient of the
log-[likelihood](https://en.wikipedia.org/wiki/Likelihood_function) for
𝒫<sub>*θ*</sub>, *i*(*θ*) is the [expected
information](https://en.wikipedia.org/wiki/Fisher_information) matrix,
and *b*(*θ*) is the *O*(*n*<sup>−1</sup>) term in the expansion of the
[bias](https://en.wikipedia.org/wiki/Bias_of_an_estimator) of the
[maximum likelihood
estimator](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation)
of *θ* (see, for example, Cox and Snell 1968).

The above iteration defines a *quasi* Fisher scoring estimation
procedure in general, and reduces to exact [Fisher
scoring](https://en.wikipedia.org/wiki/Scoring_algorithm) for maximum
likelihood estimation when *c*(*θ*<sup>(*k*)</sup>)=0<sub>*p*</sub>,
where 0<sub>*p*</sub> is a vector of *p* zeros.

For either *c*(*θ*)=*I*<sub>*p*</sub> or
*c*(*θ*)={*i*(*θ*)}<sup>−1</sup>*j*(*θ*), where *I*<sub>*p*</sub> is the
*p* × *p* identity matrix and *j*(*θ*) is the [observed
information](https://en.wikipedia.org/wiki/Observed_information),
*θ*<sup>(∞)</sup> (if it exists) is a reduced-bias estimate, in the
sense that it corresponds an estimator with bias of smaller asymptotic
order than that of the maximum likelihood estimator (see, Firth 1993;
Kosmidis and Firth 2010). The **brglm** estimates correspond to
*c*(*θ*)=*I*<sub>*p*</sub>.

The asymptotic distribution of the reduced-bias estimators is the same
to that of the maximum likelihood estimator (see, Firth 1993 for
details). So, the reduced-bias estimates can be readily used to
calculate *z*-statistics.

------------------------------------------------------------------------

Implementation using **enrichwith**
===================================

For implementing the iteration for bias reduction, we need functions
that can compute the gradient of the log-likelihood, the observed and
expected information matrix, and *b*(*θ*) at arbitrary values of *θ*.

The **enrichwith** package can produce those functions for `glm` objects
through the `auxiliary_functions` enrichment option (to see all
available enrichment options for `glm` objects run
`get_enrichment_options(modML)`.

    enriched_modML <- enrich(modML, with = "auxiliary functions")

Let's extract those functions from the `enriched_modML` object.

    # Extract the ingredients for the quasi Fisher scoring iteration from the enriched glm object
    gradient <- enriched_modML$auxiliary_functions$score # gradient of the log-likelihood
    information <- enriched_modML$auxiliary_functions$information # information matrix
    bias <- enriched_modML$auxiliary_functions$bias # first-order bias

For the more technically minded, note here that the above functions are
specific to `modML` in the sense that they look into a special
environment for necessary objects like the model matrix, the model
weights, the response vector, and so on.

This stems from the way **enrichwith** has been implemented. In
particular, if `create_enrichwith_skeleton` is used, the user/developer
can directly implement enrichment options to enrich objects with
functions that directly depend on other components in the object to be
enriched.

The following code chunk uses `enriched_modML` to implement the quasi
Fisher scoring procedure for the analysis of the endometrial cancer
data. For `p <- length(theta_mle)` , the starting value for the
parameter vector is set to `theta_current <- rep(0, p)` , and the
maximum number of iterations to `maxit <- 100` . As stopping criterion,
we use the absolute change in each parameter value with tolerance
`epsilon <- 1e-06`

    # The quasi Fisher scoring iteration using c(theta) = identity
    for (k in seq.int(maxit)) {
        s_vector <- gradient(theta_current)
        i_matrix <- information(theta_current, type = "expected")
        b_vector <- bias(theta_current)
        step <- solve(i_matrix) %*% s_vector - b_vector
        theta_current <- theta_current + step
        # A stopping criterion
        if (all(abs(step) < epsilon)) {
            break
        }
    }
    (theta_e <- drop(theta_current))

    ## (Intercept)          NV          PI          EH 
    ##  1.91460348  1.65892018 -0.01520487 -1.37987837

The estimation procedure took 8 iterations to converge, and, as
expected, the estimates are numerically the same to the ones that
**brglm** returned.

    all.equal(theta_e, theta_brglm, check.attributes = FALSE, tolerance = epsilon)

    ## [1] TRUE

A set of alternative reduced-bias estimates can be obtained using
*c*(*θ*)={*i*(*θ*)}<sup>−1</sup>*j*(*θ*). Starting again at
`theta_current <- rep(0, p)`

    # The quasi Fisher scoring iteration using c(theta) = solve(i(theta)) %*% j(theta)
    for (k in seq.int(maxit)) {
        s_vector <- gradient(theta_current)
        i_matrix <- information(theta_current, type = "expected")
        j_matrix <- information(theta_current, type = "observed")
        b_vector <- bias(theta_current)
        step <- solve(i_matrix) %*% (s_vector - j_matrix %*% b_vector)
        theta_current <- theta_current + step
        # A stopping criterion
        if (all(abs(step) < epsilon)) {
            break
        }
    }
    (theta_o <- drop(theta_current))

    ## (Intercept)          NV          PI          EH 
    ##  1.89707065  1.72815655 -0.01471219 -1.37254188

The estimation procedure took 9 iterations to converge.

The maximum likelihood estimates and the estimates from the two variants
of the bias reduction method are

    round(data.frame(theta_mle, theta_e, theta_o), 3)

    ##             theta_mle theta_e theta_o
    ## (Intercept)     2.181   1.915   1.897
    ## NV              5.805   1.659   1.728
    ## PI             -0.019  -0.015  -0.015
    ## EH             -1.526  -1.380  -1.373

Note that the reduced-bias estimates have shrunk towards zero. This is
typical for reduced-bias estimation in binomial-response generalized
linear models (see, for example, Cordeiro and McCullagh 1991, Section 8;
Kosmidis 2007, Section 5.2; Kosmidis 2014 for shrinkage in cumulative
link models).

The corresponding *z*-statistics are

    se_theta_mle <- sqrt(diag(solve(information(theta_mle, type = "expected"))))
    se_theta_e <- sqrt(diag(solve(information(theta_e, type = "expected"))))
    se_theta_o <- sqrt(diag(solve(information(theta_o, type = "expected"))))
    round(data.frame(z_mle = theta_mle/se_theta_mle,
                     z_br_e = theta_e/se_theta_e,
                     z_br_o = theta_o/se_theta_o), 3)

    ##              z_mle z_br_e z_br_o
    ## (Intercept)  2.544  2.427  2.407
    ## NV           0.009  2.220  2.215
    ## PI          -0.799 -0.728 -0.701
    ## EH          -3.523 -3.422 -3.411

The two variants for bias reduction result in slightly different
reduced-bias estimates and *z*-statistics, though the *z*-statistics
from both variants provide some evidence for the existence of an effect
for `NV`.

------------------------------------------------------------------------

Notes
=====

A general family of bias reduction methods is described in Kosmidis and
Firth (2009).

The quasi Fisher scoring iteration that has been described here is at
the core of the [**brglm2**](https://github.com/ikosmidis/brglm2) R
package, which provides various bias reduction methods for generalized
linear models.

------------------------------------------------------------------------

References
==========

Agresti, A. 2015. *Foundations of Linear and Generalized Linear Models*.
Wiley Series in Probability and Statistics. Wiley.

Cordeiro, G. M., and P. McCullagh. 1991. “Bias Correction in Generalized
Linear Models.” *Rssb* 53 (3): 629–43.

Cox, D. R., and E. J. Snell. 1968. “A General Definition of Residuals
(with Discussion).” *Journal of the Royal Statistical Society, Series B:
Methodological* 30: 248–75.

Firth, D. 1993. “Bias Reduction of Maximum Likelihood Estimates.”
*Biometrika* 80 (1): 27–38.

Heinze, G., and M. Schemper. 2002. “A Solution to the Problem of
Separation in Logistic Regression.” *Statistics in Medicine* 21:
2409–19.

Kosmidis, I. 2007. “Bias Reduction in Exponential Family Nonlinear
Models.” PhD thesis, Department of Statistics, University of Warwick.
<http://www.ucl.ac.uk/~ucakiko/files/ikosmidis_thesis.pdf>.

———. 2014. “Improved Estimation in Cumulative Link Models.” *Journal of
the Royal Statistical Society, Series B: Methodological* 76 (1): 169–96.
doi:[10.1111/rssb.12025](https://doi.org/10.1111/rssb.12025).

Kosmidis, I., and D. Firth. 2009. “Bias Reduction in Exponential Family
Nonlinear Models.” *Biometrika* 96 (4): 793–804.
doi:[10.1093/biomet/asp055](https://doi.org/10.1093/biomet/asp055).

———. 2010. “A Generic Algorithm for Reducing Bias in Parametric
Estimation.” *Electronic Journal of Statistics* 4: 1097–1112.
doi:[10.1214/10-EJS579](https://doi.org/10.1214/10-EJS579).

———. 2011. “Multinomial Logit Bias Reduction via the Poisson Log-Linear
Model.” *Biometrika* 98 (3): 755–59.
