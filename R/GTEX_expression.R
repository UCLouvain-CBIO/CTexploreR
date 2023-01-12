#' Gene expression in normal tissues (GTEx)
#'
#' @description Plots an expression heatmap of genes in normal
#' tissues (GTEx database)
#'
#' @param genes Genes selected (all CT genes by default)
#'
#' @param units Expression values unit.
#' Can be "TPM" (default) or "log_TPM" (log(TPM + 1))
#'
#' @param return If return is set to TRUE, values are returned instead of the
#' heatmap (FALSE by default)
#'
#' @param database GTEX_data
#'
#' @return A heatmap of selected genes expression in normal tissues.
#' If return = TRUE, expression values are returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' GTEX_expression(units = "log_TPM")
#' GTEX_expression(genes = c("MAGEA1", "MAGEA3"), units = "log_TPM")
GTEX_expression <- function(genes = NULL, units = "TPM", return = FALSE) {

  database <- GTEX_data
  if (is.null(genes)) genes <- CT_genes$external_gene_name
  genes <- check_gene_names(genes = genes, database = database)
  database <- database[rowData(database)$external_gene_name %in% genes, ]
  mat <- assay(database)
  rownames(mat) <- rowData(database)$external_gene_name
  name <- "TPM"
  if (units == "log_TPM") {
    mat <- log1p(mat)
    name <- "log_TPM"
  }

  if (dim(mat)[1] > 100) fontsize <- 4
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) fontsize <- 5
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) fontsize <- 6
  if (dim(mat)[1] <= 20) fontsize <- 8

  h <- suppressMessages(
    Heatmap(mat,
            name = name,
            column_title = "Gene Expression in normal tissues (GTEx)",
            col = colorRamp2(seq(0, max(mat), length = 11),
                             c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                               "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                               "#F46D43", "#D53E4F", "#9E0142")),
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            row_names_gp = gpar(fontsize = fontsize),
            column_names_gp = gpar(fontsize = 10),
            clustering_method_rows = "ward.D"))

  if (return == FALSE) {
    print(h)
  } else {
    mat
  }

}

