#' Gene expression in human embryonic stem cells
#'
#' @description
#'
#' Plots a heatmap of genes expression in human embryonic stem cells,
#' using scRNAseq data downloaded from Encode database.
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'     
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#' 
#' @param units `character(1)` with expression values unit. Can be
#' `"TPM"` (default) or `"log_TPM"` (log(TPM + 1)).
#'
#' @param scale_lims `vector of length 2` setting the lower and upper limits
#' of the heatmap colorbar. By default, the lower limit is 0, and the upper
#' limit corresponds to the third quartile of the logcounts values.
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#' SingleCellExperiment instead of the heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected CT genes expression in single cells from human
#' embryonic stem cells. If `values_only = TRUE`, a SummarizedExperiment 
#' is returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom MatrixGenerics rowMaxs
#' 
#' @examples
#' \dontrun{
#' hESC_expression(include_CTP = FALSE, units = "log_TPM",
#'                   values_only = FALSE)
#' }
hESC_expression <- function(
    genes = NULL, include_CTP = FALSE, 
    units = c("TPM"), 
    scale_lims = NULL, values_only = FALSE) {
  
    suppressMessages({
        database <- CTdata::hESC_data()
    })
    
    database <- subset_database(genes, database, include_CTP)
    
    mat <- as.matrix(SummarizedExperiment::assay(database))
    if (units == "log_TPM") mat <- log1p(mat)
    rownames(mat) <- rowData(database)$external_gene_name

    df_col <- data.frame(sample = database$sample)
    rownames(df_col) <- colnames(database)

    column_ha_sample <- HeatmapAnnotation(
      sample = df_col$sample,
        border = TRUE,
        col = list(sample = hESC_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    fontsize <- set_fontsize(mat)

    if (is.null(scale_lims)) scale_lims <- c(0, quantile(rowMaxs(mat), 0.75))

    h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
        name = units,
        column_title = "Expression in human embryonic stem cells",
        column_split = df_col$type,
        show_column_names = FALSE,
        show_column_dend = FALSE,
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        col = colorRamp2(seq(scale_lims[1], scale_lims[2], length = 11),
                         legend_colors),
        top_annotation = column_ha_sample,
        heatmap_legend_param = legends_param)

    ifelse(values_only, return(database), return(h))
}

