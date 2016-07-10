#' Build a enrichwith skeleton
#'
#'
#' Build an enrichwith skeleton file to help with the development of
#' enrichment functions and methods for objects of a specific class
#'
#' @param class the class of the objects to be enriched
#' @param option a character vector with components the enrichment
#'     options
#' @param description a character vector of length
#'     \code{length(options)} with components the description of the
#'     enrichment options
#' @param component a list of as many character vectors as
#'     \code{length(option)}, specifying the names of the components
#'     that each option will add to the object after enrichment
#' @param path the path where the skeleton file will be created
#' @param filename the name of the skeleton file
#' @param attempt_rename attempt to rename syntactically incorrect
#'     component names? Default is \code{TRUE}
#'
#' @return A file with the necessary functions to use
#'     \code{enrichwith} infrastructure. The skeleton consists of the
#'     following functions
#' \itemize{
#'
#' \item One function \code{compute_component.class} for each
#' component name from \code{unique(unlist(component))}. This function
#' takes as input the object to be enriched and returns as output the
#' component to be added to the object
#'
#' \item The \code{get_enrichment_options.class} function, that takes as
#' input the object to be enriched and returns and an enrichment
#' option, and returns the names of the components that will be
#' appended to the object for this option. This function can also be
#' used to list the available options with descriptions
#'
#' \item The \code{enrich.class} function
#'
#' }

build_enrichwith_skeleton <- function(class,
                                      option,
                                      description,
                                      component,
                                      path,
                                      filename = paste0(class, "_options.R"),
                                      attempt_rename = TRUE) {

    con <- file(paste(path, filename, sep = "/"), open = "w")

    o_length <- length(option)
    d_length <- length(description)

    message("* Checking supplied class for validity.")
    class_input <- as.character(class)
    ## Just making sure that a syntactically correct name is used for the class
    class_valid <- make.names(class, unique = TRUE, allow_ = FALSE)

    message("* Checking supplied class for validity.")
    option <- as.character(option)

    message("* Checking supplied descriptions for validity.")
    description <- as.character(description)
    if (d_length != o_length) {
        stop("ength(description) is not equal to length(option)")
    }

    message("* Checking supplied component for validity.")
    if (length(component) != o_length) {
        stop("length(component) is not equal to length(option)")
    }
    ## So even if there is only one component, component will become a list
    component <- lapply(component, function(comp) {
        comp_input <- as.character(comp)
        comp_valid <- make.names(comp, unique = TRUE, allow_ = FALSE)
        valid <- comp_input == comp_valid
        if (!all(valid)) {
            if (attempt_rename) {
                warning(gettextf("not syntactically valid component names: %s were renamed to %s",
                                 paste0(paste0("'", comp_input[!valid], "'"), collapse = ", "),
                                 paste0(paste0("'", comp_valid[!valid], collapse = ", "))))
            }
            else {
                close(con)
                stop(gettextf("not syntactically valid component names: %s",
                              paste0(paste0("'", comp_input[!valid], "'"), collapse = ", ")))
            }
        }
        comp_valid
    })
    dat <- list(class = class_input,
                option = gsub("\"","'", deparse(option, width.cutoff = 500)),
                description = gsub("\"","'", deparse(description, width.cutoff = 500)),
                component = gsub("\"","'", deparse(component, width.cutoff = 500)))

    ## enrich function
    message("* Setting up enrich.", class_input, " function")
    ## Get _options template
    template_path <- system.file("templates", "template_enrich.R", package = "enrichwith",  mustWork = TRUE)
    ## Read template and replace according to class/options/description/component
    template_out <- whisker::whisker.render(readLines(template_path), data = dat)
    ## Write into _options file
    writeLines(template_out, con = con)

    ## get_enrichment_options fuction
    message("* Setting up get_enrichment_options.", class_input, " function")
    ## Get _options template
    template_path <- system.file("templates", "template_get_enrichment_options.R", package = "enrichwith",  mustWork = TRUE)
    ## Read template and replace according to class/options/description/component
    template_out <- whisker::whisker.render(readLines(template_path), data = dat)
    ## Write into _options file
    writeLines(template_out, con = con)

    ## compute_component functions
    message("* Setting up compute_component.", class_input, " function")
    components <- unique(unlist(component))
    template_path <- system.file("templates", "template_compute_component.R", package = "enrichwith",  mustWork = TRUE)
    for (j in seq.int(length(components))) {
        template_out <-  whisker::whisker.render(readLines(template_path),
                                                 data = list(component = components[j],
                                                             class = class_input))
        writeLines(template_out, con = con)
    }

    close(con)
}


