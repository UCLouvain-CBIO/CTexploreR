#' Keeps genes present in a specified SE
#'
#' @description Checks if specified genes are present in the rowData's column 
#' `external_gene_name` of a SummarizedExperiment, and returns the genes that 
#' are in it. 
#'
#' @param genes Genes 

#' @param database `SE` The summarizeExperiment must contain a column 
#' external_gene_name in it rowData.
#'
#' @return genes present in the specified database
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @examples
#' check_gene_names(genes = c("MAGEA1", "MAGEA3", "XXX"), database = GTEX_data)
check_gene_names <- function(genes, database) {
  if (!all(genes %in% rowData(database)$external_gene_name)) {
    message("Check gene name(s)!\n")
    message(paste0(
      genes[!genes %in% rowData(database)$external_gene_name],
      " is not in the database.\n"))
  }
  genes <- genes[genes %in% rowData(database)$external_gene_name]
}
