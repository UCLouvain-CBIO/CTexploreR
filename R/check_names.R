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
#' check_names(variable = c("Ovarian", "leukemia", "wrong_name"),
#'             valid_vector = c("ovarian", "leukemia"))
check_names <- function(variable, valid_vector) {
  if (!all(variable %in% valid_vector)) {
      message(paste(variable[!variable %in% valid_vector],
                    "is not valid: check spelling!\n"))

      if (length(valid_vector) < 20) {
          message("Should be one of: ",
                  paste(valid_vector, collapse = " "),
                  "\n")
    }
      variable <- variable[variable %in% valid_vector]
      if (length(variable) > 0 & length(variable) < 20) {
          message("Kept only: ", paste(variable, collapse = " "))
    }
  }
  variable
}
