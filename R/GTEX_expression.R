#' Gene expression in normal tissues (GTEx)
#'
#' @description Plots an expression heatmap of genes in normal
#' tissues (GTEx database)
#'
#' @param database GTEX_data
#'
#' @param genes Genes selected (all CT genes by default)
#'
#' @param units Expression values unit.
#' Can be "TPM" (default) or "log_TPM" (log(TPM + 1))
#'
#' @return A heatmap of selected genes in normal tissues.
#' Expression values are invisibly returned.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' GTEX_expression(database = GTEX_data)
#' GTEX_expression(database = GTEX_data, genes = c("MAGEA1", "MAGEA3"),
#' units = "log_TPM")
GTEX_expression <- function(database, genes = NULL, units = "TPM") {

  if (missing(database)) {
    stop("Database must be specified!")
  }

  if (!missing(database)) {
    mat <- assay(database)
    rownames(mat) <- rowData(database)$external_gene_name
  }

  if (is.null(genes)) {
    mat <- mat[CT_genes$external_gene_name,]
  }

  if (!is.null(genes)) {
    if (!all(genes %in% rownames(mat))) { ## Warning !
      message("Check gene name(s)!\n")
      message(paste0(
        genes[!genes %in% rownames(mat)],
        " is not in the database.\n"))
      genes <- genes[genes %in% rownames(mat)]
    }
    mat <- mat[genes, ]
  }

  name <- "TPM"
  if (units == "log_TPM") {
    mat <- log1p(mat)
    name <- "log_TPM"
  }

  if (dim(mat)[1] > 100){ fontsize <- 4 }
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100){ fontsize <- 5 }
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50){ fontsize <- 6 }
  if (dim(mat)[1] <= 20) { fontsize <- 8 }

  h <- suppressMessages(Heatmap(mat,
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
  print(h)
  invisible(mat)
}

