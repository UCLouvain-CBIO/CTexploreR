#' Gene expression in different human cell types
#'
#' @description
#'
#' Plots a heatmap of genes expression in the different human cell types based
#' on scRNAseq data obtained from the Human Protein Atlas
#' (https://www.proteinatlas.org)
#'
#' @param genes `character` naming the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'     
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).  
#'
#' @param units `character(1)` with expression values unit.  Can be
#'     `"TPM"`, `"log_TPM"` (log(TPM + 1)) or `"scaled"` (scaled TPM
#'     values). Default is `"scaled"`.
#'
#' @param scale_lims `vector of length 2` setting the lower and upper limits
#' of the heatmap colorbar.
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#'     SummarizedExperiment instead of the heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected CT genes expression in
#' different human cell types.
#' If `values_only = TRUE`, a SummarizedExperiment instead of the heatmap
#'     is returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom MatrixGenerics rowMaxs
#'
#' @examples
#'
#' HPA_cell_type_expression(
#'     genes = NULL, units = "scaled", scale_lims = NULL,
#'     values_only = FALSE)
#' HPA_cell_type_expression(
#'     genes = c("MAGEA1", "MAGEA3", "MAGEA4"),
#'     units = "TPM", scale_lims = c(0, 50),
#'     values_only = FALSE)
HPA_cell_type_expression <- function(genes = NULL, 
                                     units = c("scaled", "TPM", "log_TPM"),
                                     include_CTP = FALSE, 
                                     scale_lims = NULL, values_only = FALSE) {
    suppressMessages({
        database <- CTdata::scRNAseq_HPA()
    })
  
    units <- match.arg(units)
  
    database <- subset_database(genes, database, include_CTP)

    ## Use gene names instead of ENSEMBL IDs
    mat <- SummarizedExperiment::assay(database)
    rownames(mat) <- rowData(database)$external_gene_name

    df_col <- data.frame(group = database$group)
    rownames(df_col) <- database$Cell_type
    not_somatic_group <- c("Germ_cells", "Trophoblast_cells")
    somatic_groups <- unique(database$group[!database$group %in%
        not_somatic_group])
    df_col$group <- factor(df_col$group,
        levels = c(not_somatic_group, somatic_groups))

    column_ha_group <- HeatmapAnnotation(
        group = df_col$group,
        border = TRUE,
        col = list(group = cell_type_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    fontsize <- set_fontsize(mat)

    if (units == "log_TPM") mat <- log1p(mat)
    if (units == "scaled") mat <- na.omit(t(scale(t(mat))))
    if (is.null(scale_lims)) scale_lims <- c(0, rowMaxs(mat))

    h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
        name = units,
        column_title = "Expression in human cell types",
        column_split = df_col$group,
        show_column_names = TRUE,
        show_column_dend = FALSE,
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = 8),
        col = colorRamp2(seq(scale_lims[1], scale_lims[2], length = 11),
                         legend_colors),
        top_annotation = column_ha_group,
        heatmap_legend_param = legends_param)

    ifelse(values_only, return(database), return(h))
}
