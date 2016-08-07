# enrichwith

Methods to enrich list-like R objects with extra components

## Installation
Get the development version from github with

``` r
# install.packages("devtools")
devtools::install_github("ikosmidis/enrichwith")
```

## Aim and objective

Suppose you developed a piece of statistical methodology that relies
on a component that a list-like R object `object_x` of a certain class
could potentially have but doesn't.

The aim of **enrichwith** is to:

* produce customisable source code templates for the structured implementation of methods to compute new components
* provide generic methods for the easy enrichment of the object with those components

Specifically, once the methods for the components have been
implemented, `object_x` can be enriched with the components
corresponding to the option `enrichment_option` through the following simple call.

``` r
enrich(object_x, with = enrichment_option)
```

The main objectives of **enrichwith** is to allow users and developers to directly use the enrichment options that other developers have provided, minimising the need of adapting source code of others.

## Example
Objects of class `link-glm` have as components functions to compute the link function (`linkfun`), the inverse link function (`linkinv`), and the 1st derivative of the link function (`mu.eta`).

**enrichwith** comes with a built-in template with the methods for the enrichment of `link-glm` objects with the 2nd and 3rd derivatives of the inverse link function.

The `get_enrichment_options` method can be used to check what enrichment options are available for objects of class `link-glm`.
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

According to the result of `get_enrichment_options`, the object `standard_link` can be enriched with the 2nd and 3rd derivative of the inverse
link function through the option "inverse link derivatives".
``` r
enriched_link <- enrich(standard_link, with = "inverse link derivatives")
class(enriched_link)
# [1] "link-glm"
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
`enriched_link` is now an "enriched" `link-glm` object, which, as per the enrichment options above, has the extra components `d2mu.deta` and `d3mu.deta`, for the calculation of 2nd and 3rd derivatives of the inverse link function with respect to `eta`, respectively.

## Implementation of enrichment options
The task of implementing the enrichment options is streamlined into 3 steps:

1. Use `create_enrichwith_skeleton` to produce an enrichwith template.
2. Add the specific code that calculates the components by editing the
   `compute_*` functions.
3. Finalise the documentation and/or include more examples.

The first step results in a template that includes all necessary functions to carry out the enrichment. The second step is where the user edits that template and adds code to calculate the components that the object will be enriched with. Specifically, each `compute_*` function takes as input the object to be enriched and returns the corresponding new component to be added to the object.

Everything else (for example, mapping between the enrichment options and the components that the enriched object will have, checks that an enrichment option exists, listing enrichment options, enriching the object, and so on) is taken care of by the methods in **enrichwith**.

Developers can either put their **enrichwith** templates in their packages or are welcome to contribute their template to **enrichwith**, particularly if that extends core R objects.

