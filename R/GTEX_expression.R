#' Gene expression in normal tissues (GTEx)
#'
#' @description Plots an expression heatmap of genes in normal tissues
#'     (GTEx database).
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#' 
#' @param units `character(1)` with expression values unit.  Can be
#'     `"TPM"` (default) or `"log_TPM"` (log(TPM + 1)).
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#'     expression values in all samples instead of the
#'     heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected genes expression in normal tissues.
#'     If `values_only = TRUE`, expression values are returned instead.
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
GTEX_expression <- function(genes = NULL, units = c("TPM", "log_TPM"), 
                            include_CTP = FALSE, values_only = FALSE) {
    
  suppressMessages({
    database <- CTdata::GTEX_data()
    })
  
    units <- match.arg(units)

    database <- subset_database(genes, database, include_CTP)

    mat <- assay(database)
    rownames(mat) <- rowData(database)$external_gene_name
    
    if (units == "log_TPM") mat <- log1p(mat)

    fontsize <- set_fontsize(mat)

    h <- Heatmap(mat,
        name = units,
        column_title = "Gene Expression in normal tissues (GTEx)",
        col = colorRamp2(seq(0, max(mat), length = 11), legend_colors),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = 10),
        clustering_method_rows = "ward.D")

    ifelse(values_only, return(mat), return(h))
}
