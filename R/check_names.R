#' Check spelling of entered variables
#'
#' @description Checks the spelling of a vector of entered variable(s) comparing
#' it to a vector of valid names, and removes the ones that are absent from the
#' vector of valid names.
#'
#' @param variable `vector` Names of variables to check
#'
#' @param valid_vector `vector` Valid variable names
#'
#' @return A vector with valid variables
check_names <- function(variable, valid_vector) {

  if (!all(variable %in% valid_vector)) {
    message(paste(variable[!variable %in% valid_vector], "is not valid: check spelling!\n"))
    if (length(valid_vector) < 20){
    message("Should be one of :")
    message(paste0(valid_vector, " "))
    }
    variable <- variable[variable %in% valid_vector]
    if (length(variable) > 0 & length(variable) < 20){
      message("\nKept only:")
      message(paste0(variable, " "))
    }
  }
  variable
}


