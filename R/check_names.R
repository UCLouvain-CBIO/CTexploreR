#' Check spelling of entered variables
#'
#' @description
#'
#' Checks the spelling of a vector of entered variable(s) comparing it
#' to a vector of valid names, and removes the ones that are absent
#' from the vector of valid names.
#'
#' @param variable `character()` containing the names of variables to
#'     check.
#'
#' @param valid_vector `character()` with valid variable names.
#'
#' @return A character with valid variables.
#'
#' @examples
#' CTexploreR:::check_names(variable = c("Ovarian", "leukemia", "wrong_name"),
#'                          valid_vector = c("ovarian", "leukemia"))
check_names <- function(variable, valid_vector) {
    in_valid <- variable %in% valid_vector
    if (!all(in_valid)) {
        msg <- paste0(sum(!in_valid), " out of ",
                     length(in_valid), " names invalid: ",
                     paste(variable[!in_valid], collapse = ", "),
                     ".")
        warning(paste(strwrap(msg), collapse = "\n"),
                "\nSee the manual page for valid types.",
                call. = FALSE)
    }
    variable[in_valid]
}
