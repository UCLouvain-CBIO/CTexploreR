#' Gene expression in CCLE Tumors
#'
#' @description
#'
#' Plots an expression heatmap of genes in CCLE tumor cell lines.
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT genes.
#'
#' @param type `character` describing the tumor cell line(s) type to
#'     be plotted. Allowed cell lines are "Ovarian", "Leukemia",
#'     "Colorectal", "Skin", "Lung", "Bladder", "Kidney", "Breast",
#'     "Pancreatic", "Myeloma", "Brain", "Sarcoma", "Lymphoma",
#'     "Bone", "Neuroblastoma", "Gastric", "Uterine", "Head_and_Neck",
#'     "Bile_Duct" and "Esophageal".
#'
#' @param type `character()` describing the tumor cell line(s) type to
#'     be plotted. Allowed cell lines are "Ovarian", "Leukemia",
#'     "Colorectal", "Skin", "Lung", "Bladder", "Kidney", "Breast",
#'     "Pancreatic", "Myeloma", "Brain", "Sarcoma", "Lymphoma",
#'     "Bone", "Neuroblastoma", "Gastric", "Uterine", "Head_and_Neck",
#'     "Bile_Duct" and "Esophageal".
#'
#' @param units `character(1)` with expression values unit. Can be
#'     "TPM" (default) or "log_TPM" (log(TPM + 1))
#'
#' @param return `logical(1)`. If `TRUE`, values are returned instead
#'     of the heatmap (`FALSE` by default).
#'
#' @return A heatmap of selected genes in CCLE cell lines from
#'     specified type.  If `return` is `TRUE`, expression values are
#'     returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @import CTdata
#'
#' @examples
#' CCLE_expression(
#'     genes = c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEA6", "MAGEA10"),
#'     type = c("Skin", "Lung"), units = "log_TPM")
CCLE_expression <- function(genes = NULL, type = NULL, units = "TPM",
                            return = FALSE) {
    suppressMessages({
        database <- CTdata::CCLE_data()
        CT_genes <- CTdata::CT_genes()
    })

    database$type <- tolower(database$type)
    valid_tumor_types <- unique(database$type)
    type <- check_names(variable = tolower(type),
                        valid_vector = valid_tumor_types)
    stopifnot("No valid tumor type entered" = length(type) > 0)
    database <- database[, database$type %in% type]

    database <- subset_database(genes, database)

    mat <- assay(database)
    rownames(mat) <- rowData(database)$external_gene_name

    name <- "TPM"
    if (units == "log_TPM") {
        mat <- log1p(mat)
        name <- "log_TPM"
    }

    df_col <- data.frame(
        "cell_line" = colData(database)$cell_line_name,
        "type" = as.factor(colData(database)$type))
    rownames(df_col) <- rownames(colData(database))
    df_col <- df_col[order(df_col$type), ]

    column_ha_type <- HeatmapAnnotation(
        type = df_col$type,
        border = TRUE,
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param,
        col = list(type = CCLE_colors))

    fontsize <- setFontSize(mat)
    
    if (length(type) <= 5) label_fontsize <- 6
    if (length(type) > 5 & length(type) < 10) label_fontsize <- 4
    if (length(type) >= 10 | is.null(type)) label_fontsize <- 0

    h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
        name = name,
        column_title = "Gene Expression in tumor cell lines (CCLE)",
        column_split = factor(df_col$type),
        col = colorRamp2(seq(0, max(mat), length = 11),legend_colors),
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        show_column_names = FALSE,
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        heatmap_legend_param = legends_param,
        top_annotation = c(column_ha_type))

    if (return) {
        return(mat)
    }
    return(h)
}