# enrichwith 0.3

## Bug fixes
* Added conditional use of `SuppDists`

## Other improvements, updates and additions
* Fixed typos and spelling mistakes in documentation
* Added new tests for d/r/p/qmodel

# enrichwith 0.3

## New functionality
* Added support for `wedderburn` family from the `gnm` R package.

## Other improvements, updates and additions
* Improved documentation
* Improved tests and added some new
 
# enrichwith 0.2

## Bug fixes
* Fixed bug in the computation of beta regression score contributions from `betareg` obejcts
* Fixed bug with the attributes of information matrices for generalized linear models with fixed dispersion
* Fixed the row names in the information matrices for `glm` and `betareg` objects

## New functionality
* `get_information_function` for `betareg` objects returns a function that can compute observed information matrices (by setting `type = "observed"`)
* `get_information_function` now returns a function that can compute the Cholesky decomposition of the information matrix (by setting `CHOL = TRUE`)
* `simulate` methods (coming form `get_simulate_function`) now check for correct length of coefficient vector if this is specified

## Other improvements, updates and additions
* Fixed various problems with documentation
* Improved tests and added some new
* Changed maintainer email to <ioannis.kosmidis@warwick.ac.uk>

# enrichwith 0.1

## Bug fixes

## New functionality
* `enrich.family` can return 4th derivatives of the a function and other characteristics of the exponential family
* `enrich` method for `betareg` objects (enrichment options with scores, information, bias and simulate)

## Other improvements, updates and additions
* Documentation updates
* Added vignette for `enrich.family` including a description of the exponential family

# enrichwith 0.05

## Bug fixes
* Fixed bug when passing a vector of options in the `with` argument of `enrich` methods
* Fixed a bug that would return wrong results if the `bias` function was passed a named vectors

## New functionality

## Other improvements, updates and additions
* Added d,p,q and enriched_glm in the vignette
* Documentation updates

# enrichwith 0.04

## Bug fixes
* p,q,d model: if data is missing then the model frame is used
* Fixed error message in qmodel

## New functionality
* Added `get_dmodel_function`, `get_pmodel_function`, `get_qmodel_function`
* `enriched_glm` can now be used to fit GLMs and get objects that are
  enriched with auxiliary functions and other components when compared
  to their `glm` counterparts

## Other improvements, updates and additions
* Fixed typos in `?enrich.family`
* Added documentation for `get_dmodel_function.glm`,
  `get_pmodel_function.glm`, `get_qmodel_function.glm`

# enrichwith 0.03

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

# enrichwith 0.02

* Fixed an issue with the example in ?enrich.glm
* Minor codebase improvements
* Harmonised the output of the functions in the auxiliary_functions component of enriched `glm` objects
* More detailed descriptions for the enrichment options for `glm` objects
* Included a vignette on how bias-reduction for GLMs can be implemented using enriched `glm` objects

# enrichwith 0.01

* First release



