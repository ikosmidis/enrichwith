# enrichwith 0.3

## Bug fixes
* Took care of aliasing in bias calculations
* Names of score and bias components do not get a prefix anymore
* Fixed a bug on argument of information function (part of the "auxiliary functions" enrichment option)
* Fixed a bug when testing the simulation methods

## New functionality
* Provide enrich capabilities for `lm` objects
* Score functions (?get_score_function) now accept `contributions = TRUE` for getting score contributions
* Provide the get_*_function convenience methods (e.g. get_score_function, get_information_function, get_bias_function)
* Provide the dmodel, pmodel, qmodel auxiliary functions for glm objects, to compute d, p, q based on a supplied data frame

## Other improvements, updates and additions
* Various codebase improvements
* Various documentation improvements
* New vignette for enriching `glm` objects
* Updated README file

# enrichwith 0.2

* Fixed an issue with the example in ?enrich.glm
* Minor codebase improvements
* Harmonised the output of the functions in the auxiliary_functions component of enriched `glm` objects
* More detailed descriptions for the enrichment options for `glm` objects
* Included a vignette on how bias-reduction for GLMs can be implemented using enriched `glm` objects

# enrichwith 0.1

* First release



