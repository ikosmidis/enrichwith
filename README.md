# enrichwith

## Goal
**enrichwith** provides a framework for enriching R objects with extra components.

## Installation

Get the development version from github with

``` r
# install.packages("devtools")
devtools::install_github("ikosmidis/enrichwith")
```

## Examples

Objects of class `link-glm` have as components functions that can be used to evaluate the link function (`linkfun`), the inverse link function (`linkinv`), and the 1st derivative of the link function (`mu.eta`).

Here we use **enrichwith** to enrich `link-glm` objects with the 2nd and 3rd derivatives of the inverse link function.

First, let's use the `get_enrichment_options` method to check what enrichment options are available for objects of class `link-glm`.
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

Below we use the "inverse link derivatives" option to enrich the
object `standard_link` with the 2nd and 3rd derivative of the inverse
link function:
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
`enriched_link` is now an enriched `link-glm` object, which, as per the enrichment options above, has the extra components `d2mu.deta` and `d3mu.deta`, for the calculation of 2nd and 3rd derivatives of the inverse link function with respect to `eta`, respectively.

## Details
The **enrichwith** package is an attempt to streamline the process of enriching objects with new components.

For example, suppose you developed a piece of statistical methodology that relies on a property (e.g. 2nd derivatives of the inverse link function) that an object of certain class (e.g. `link-glm`) could have but doesn't.

**enrichwith** provides methods that enable the implementation of the new property, and the enrichment of an object with the corresponding component in a *principled* way. The resulting object

* has the same components as the original object plus the extra compoments from the enrichment;
* preserves the class of the original object, so that all methods that were associated with the original object, are also associated with the new object.

The task of implementing the enrichment options is streamlined into 3 steps:

1. Use `create_enrichwith_skeleton` to produce an enrichwith template.
2. Add the specific code that calculates the components by editing the
   `compute_*` functions.
3. Finalise the documentation and/or include more examples.

The first step results in a template that includes all necessary functions to carry out the enrichment. The second step is where the user edits that template and adds code to calculate the components that the object will be enriched with. Specifically, each `compute_*` function takes as input the object to be enriched and returns the corresponding new component to be added to the object.

Everything else (for example, mapping between the enrichment options and the components that the enriched object will have, checks that an enrichment option exists, listing enrichment options, enriching the object, and so on) is taken care of by the methods in **enrichwith**.
