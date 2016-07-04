#' Return a template for contruction a "_options" function
#'
#' @param class the class of the object for which options are to be listed
#'
#' @details The function prints a new function that can be used to
#'     define and check for the available options for enrichment of
#'     the object with class \code{class}. The user is expected to
#'     populate the \code{available_options} vector and the
#'     \code{*_with} list. The latter may simply be
#'     \code{as.list(available_options)} if no groups of options are
#'     needed.
#'
#' @return A function.
#'
#' @examples
#' \dontrun{
#' ## family_options has been written via
#' options_function("family")
#' }
options_function <- function(class) {
    cat(paste0(class, "_options <- function(what, list_options = FALSE) {"), "\n")
    cat("\t", "## List the enrichment options that you would like to make\n")
    cat("\t", "## avaiable for objects of class\n")
    cat("\t", "available_options <- c(!!COMPLETEME!!))", "\n")
    cat("\t", "## Provide the descriptions of the enrichment options\n")
    cat("\t", "descriptions <- c(!!COMPLETEME!!))", "\n")
    cat("\t", "available_options <- c(available_options, 'all')\n")
    cat("\t", "descriptions <- c(descriptions, 'all available options')\n")
    cat("\t", "if (list_options) {\n")
    cat("\t", "cat(paste(paste(available_options, descriptions, sep = ' : '), '\\n'))\n")
    cat("\t", "return(invisible())\n")
    cat("\t", "}", "\n")
    cat("\t", "if (any(!(what %in% available_options))) {\n")
    cat("\t", "stop(gettextf('one of the options %s is not implemented', paste(what, collapse = ', ')))\n")
    cat("\t", "}\n")
    cat("\t", "## List what each option in available_options corresponds to\n")
    cat("\t", "## (one vector of function names per option)l. The corresponding\n")
    cat("\t", "## functions should take as input an object of class 'class'\n")
    cat("\t", paste0(class, "_with <- list(!!COMPLETEME!!)"), "\n")
    cat("\t", paste0(class, "_with[[length(", class, "_with) + 1]] <- unique(unlist(", class, "_with))", collapse = ""), "\n")
    cat("\t", paste0("names(", class, "_with) <- available_options"), "\n")
    cat("\t", paste0("unique(unlist(", class, "_with[what]))"), "\n")
    cat("}", "\n")
}
