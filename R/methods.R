#' @export
print.enrichment_options <- function(x, ...) {
    o_length <- length(x$option)
    for (o in seq.int(o_length)) {
        cat('-------\n')
        cat('Option:', x$option[o], "\n")
        cat('Description:', x$description[o], "\n")
        mat <- data.frame(component = x$component[[o]], compute_function = x$compute_function[[o]])
        print(mat)
    }
}
