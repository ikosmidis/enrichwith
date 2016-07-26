enrichwith
========

### Details

The **enrichwith** package provides the `enrich` method to enrich various core objects with extra, relevant components. The resulting objects preserve their class, so all methods associated to them still apply.

Depending on the object, enriching it can be a tedious task. The
**enrichwith** package is an attempt to structure the task into 3 simple steps:

1. Use `create_enrichwith_skeleton` to produce an enrichwith template.
2. Complete the functions for the `compute_*` methods in the enrichwith template.
3. Finalise the documentation and/or include more examples.


The first step results in a template that includes all necessary functions to carry out the enrichment. The second step is where the user writes the code to calculate the components that the object will be enriched with. Specifically, each `compute_*` function takes as input the object to be enriched and return the corresponding new component to be used for the enrichment of the object.

Everything else (for example, mapping between the enrichment options and the components that the enriched object will have, checks that an enrichment option exists, listing enrichment options, enriching the object, and so on) is taken care of by the methods in **enrichwith**.

### Installation

Install the development version from github:

``` r
devtools::install_github("ikosmidis/enrichwith")
```

