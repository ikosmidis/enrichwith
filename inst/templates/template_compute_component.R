`compute_{{component}}.{{class}}` <- function(object, ...) {
    ## Write some code to compute the component {{component}} using
    ## the components of object and any other arguments in ...
    ## For example
    cat('component', '{{component}}', 'for objects of class', '{{class}}', "\n")
}


`compute_{{component}}` <- function(object, ...) {
    UseMethod('compute_{{component}}')
}


