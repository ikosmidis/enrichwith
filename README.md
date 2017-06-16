[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/enrichwith)](https://cran.r-project.org/package=enrichwith)

# enrichwith

Methods to enrich list-like R objects with extra components

## Rationale

When I develop a piece of statistical methodology it is not uncommon
that it depends on functionality or quantities that a list-like core R
object `object_x` of a certain class could have but doesn't. An
example is objects of class `link-glm` which only provide up to the
the first derivatives of the link function. It is also not uncommon
that the functionality or quantities that I need are hard-coded in
another package in the R ecosystem.  Prominent examples include
implementations of gradients of the log-likelihood and information
matrices for specific model classes.

In such cases, I either contact the developers and ask them to provide
a generic which I can then use, or copy some of the code and adopt it
to what I'm doing. Both options can be rather time-consuming, and
particularly the latter is rarely bug-free.

I believe that

* users and developers should have direct access to useful
  functionality and quantities in the R ecosystem, epsecially if these include implementations of complex statistical quantities

* quantities and functionality that are specific to a list-like `object_x` should be components of `object_x`

The above motivated me to develop the **enrichwith** R package that
allows `object_x` to be enriched with components corresponding to the
option `enrichment_option` through the following simple call

``` r
enrich(object_x, with = enrichment_option)
```

The call is inspired by [Donald
Knuth's](https://en.wikipedia.org/wiki/Donald_Knuth) [literare
programing](https://en.wikipedia.org/wiki/Literate_programming)
paradigm.

## Purpose and objective

The main objective of **enrichwith** is to allow users and developers
to directly use the enrichment options that other developers have
provided, through a *clean interface*, minimising the need to adopt
source code of others.

The purpose of **enrichwith** is to provide:

* *useful enrichment options* for core R objects, including objects of
  class `lm`, `glm`, `link-glm` and `family` (see, for example,
  `?enrich.glm`)

* methods for producing *customisable source code templates* for the
  structured implementation of methods to compute new components (see
  `?enrichwith` and `?create_enrichwith_skeleton`)

* generic methods for the *easy enrichment* of the object with those
  components (see, for example, `?enrich` and
  `?get_enrichment_options`)

## Vignettes

The vignettes illustrate the **enrichwith** functionality through
comprehensive, step-by-step case studies. These also include
illustrations from recent research of mine on methods for statistical
learning and inference.

## Installation
Get the development version from github with

``` r
# install.packages("devtools")
devtools::install_github("ikosmidis/enrichwith")
```

## Example

Objects of class `link-glm` have as components functions to compute
the link function (`linkfun`), the inverse link function (`linkinv`),
and the 1st derivative of the inverse link function (`mu.eta`).

**enrichwith** comes with a built-in template with the methods for the
enrichment of `link-glm` objects with the 2nd and 3rd derivatives of
the inverse link function.

The `get_enrichment_options` method can be used to check what
enrichment options are available for objects of class `link-glm`.
``` r
library("enrichwith")
standard_link <- make.link("probit")
class(standard_link)
# [1] "link-glm"
get_enrichment_options(standard_link)
# -------
# Option: d2mu.deta
# Description: 2nd derivative of the inverse link function
#   component  compute_function
# 1 d2mu.deta compute_d2mu.deta
# -------
# Option: d3mu.deta
# Description: 3rd derivative of the inverse link function
#   component  compute_function
# 1 d3mu.deta compute_d3mu.deta
# -------
# Option: inverse link derivatives
# Description: 2nd and 3rd derivative of the inverse link function
#   component  compute_function
# 1 d2mu.deta compute_d2mu.deta
# 2 d3mu.deta compute_d3mu.deta
# -------
# Option: all
# Description: all available options
#   component  compute_function
# 1 d2mu.deta compute_d2mu.deta
# 2 d3mu.deta compute_d3mu.deta
```

According to the result of `get_enrichment_options`, the object
`standard_link` can be enriched with the 2nd and 3rd derivative of the
inverse link function through the option "inverse link derivatives".
``` r
enriched_link <- enrich(standard_link, with = "inverse link derivatives")
class(enriched_link)
# [1] "enriched_link-glm" "link-glm"
cat(format(enriched_link$d2mu.deta), sep = "\n")
# function (eta)
# {
#     -eta * pmax(dnorm(eta), .Machine$double.eps)
# }
cat(format(enriched_link$d3mu.deta), sep = "\n")
# function (eta)
# {
#     (eta^2 - 1) * pmax(dnorm(eta), .Machine$double.eps)
# }
```
`enriched_link` is now an "enriched" `link-glm` object, which, as per
the enrichment options above, has the extra components `d2mu.deta` and
`d3mu.deta`, for the calculation of 2nd and 3rd derivatives of the
inverse link function with respect to `eta`, respectively.

