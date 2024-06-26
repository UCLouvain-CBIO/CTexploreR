#' Gene expression in TCGA tumors
#'
#' @description
#'
#' Plots a heatmap of genes expression in TCGA samples (peritumoral
#' and tumor samples when a specific tumor type is specified, or tumor
#' samples only when tumor option is set to "all")
#'
#' @param tumor `character` defining the TCGA tumor type. Can be one
#'     of "SKCM", "LUAD", "LUSC", "COAD", "ESCA", "BRCA", "HNSC", or
#'     "all" (default).
#'
#' @param genes `character` naming the selected genes. The default
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
#' @return A heatmap of selected CT genes expression in TCGA samples.
#'     If `values_only = TRUE`, TPM expression data is returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' \dontrun{
#' TCGA_expression(
#'     tumor = "LUAD", genes = c("MAGEA1", "MAGEA3"),
#'     units = "log_TPM")
#' }
TCGA_expression <- function(tumor = "all", genes = NULL,
                            include_CTP = FALSE,
                            units = c("TPM", "log_TPM"), values_only = FALSE) {
    suppressMessages({
        database <- CTdata::TCGA_TPM()
    })
    
    units <- match.arg(units)

    database$tumor <- sub(pattern = "TCGA-", x = database$project_id, "")
    valid_tumors <- unique(database$tumor)
    tumor <- check_names(tumor, c(valid_tumors, "all"))
    stopifnot("No valid tumor type entered" = length(tumor) > 0)
    if (!"all" %in% tumor) database <- database[, database$tumor %in% tumor]

    database$type <- "Tumor"
    database$type[database$shortLetterCode == "NT"] <- "Peritumoral"
    database <- database[, order(database$tumor, database$type)]

    database <- subset_database(genes, database, include_CTP)

    ## Peritumoral samples are displayed only when a single type of tumor asked
    if ("all" %in% tumor | length(tumor) > 1) {
        database <- database[, database$type == "Tumor"]
        split_by <- factor(database$tumor)
        column_ha_tumor <- HeatmapAnnotation(
          Tumor = database$tumor,
          border = TRUE,
          col = list(Tumor = TCGA_colors),
          annotation_name_gp = gpar(fontsize = 8),
          annotation_legend_param = legends_param)
        annot <- column_ha_tumor
    } else {
        split_by <- factor(database$type, levels = c("Peritumoral", "Tumor"))
        column_ha_type <- HeatmapAnnotation(
          Type = database$type,
          border = TRUE,
          col = list(Type = c("Tumor" = "salmon",
                              "Peritumoral" = "mediumseagreen")),
          annotation_name_gp = gpar(fontsize = 8),
          annotation_legend_param = legends_param)
        annot <- column_ha_type
    }

    ## Use gene names instead of ENSEMBL IDs
    mat <- SummarizedExperiment::assay(database)
    rownames(mat) <- rowData(database)$external_gene_name

    if (units == "log_TPM") mat <- log1p(mat)

    fontsize <- set_fontsize(mat)

    h <- Heatmap(mat[, , drop = FALSE],
        name = units,
        column_title = paste0("Expression in TCGA samples (", tumor, ")"),
        column_split = split_by,
        col = colorRamp2(seq(0, max(mat), length = 11), legend_colors),
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        cluster_rows = TRUE,
        show_column_names = FALSE,
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        heatmap_legend_param = legends_param,
        top_annotation = annot)

    ifelse(values_only, return(mat), return(h))
}
