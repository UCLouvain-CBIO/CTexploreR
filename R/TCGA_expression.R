#' Gene expression in TCGA tumors
#'
#' @description Plots a heatmap of genes expression in TCGA samples
#' (peritumoral and tumor samples when a specific tumor type is specified, or
#' tumor samples only when tumor option is set to "all")
#'
#' @param tumor TCGA tumor code. Can be one of "SKCM", "LUAD", "LUSC", "COAD",
#' "ESCA", "BRCA", "HNSC", or "all" (default).
#'
#' @param genes Genes selected (All CT genes by default)
#'
#' @param units Expression values unit.
#' Can be "TPM" (default) or "log_TPM" (log(TPM + 1))
#'
#' @param return Boolean (FALSE by default). If set to TRUE, the function will
#' return the expression values in all samples instead of the heatmap.
#'
#' @return A heatmap of selected CT genes expression in TCGA samples.
#' If return = TRUE, TPM expression data is returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' TCGA_expression(tumor = "LUAD", genes = c("MAGEA1", "MAGEA3"), units = "log_TPM")
#' TCGA_expression(tumor = "all", units = "log_TPM")
TCGA_expression <- function(tumor = "all", genes = NULL,
                            units = "TPM", return = FALSE) {
    database <- CTdata::TCGA_TPM()
    CT_genes <- CTdata::CT_genes()

    database$tumor <- sub(pattern = 'TCGA-', x = database$project_id, '')
    valid_tumors <- unique(database$tumor)
    tumor <- check_names(tumor, c(valid_tumors, "all"))
    stopifnot("No valid tumor type entered" = length(tumor) > 0)
    if (!"all" %in% tumor) database <- database[, database$tumor %in% tumor]

    database$type <- "Tumor"
    database$type[database$shortLetterCode == "NT"] <- "Peritumoral"
    database <- database[, order(database$tumor, database$type)]

    if (is.null(genes)) genes <- CT_genes$external_gene_name
    valid_gene_names <- unique(rowData(database)$external_gene_name)
    genes <- check_names(genes, valid_gene_names)
    database <- database[rowData(database)$external_gene_name %in% genes, ]

    legends_param <- list(
        labels_gp = gpar(col = "black", fontsize = 6),
        title_gp = gpar(col = "black", fontsize = 6),
        row_names_gp = gpar(fontsize = 4),
        annotation_name_side = "left")

    column_ha_type <- HeatmapAnnotation(
        Type = database$type,
        border = TRUE,
        col = list(Type = c("Tumor" = "salmon",
                            "Peritumoral" = "mediumseagreen")),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    column_ha_tumor <- HeatmapAnnotation(
        Tumor = database$tumor,
        border = TRUE,
        col = list(Tumor = TCGA_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    ## Peritumoral samples are displayed only when a single type of tumor is asked
    if ("all" %in% tumor | length(tumor) > 1) {
        split_by <- factor(database$tumor)
        annot <- column_ha_tumor
    } else {
        split_by <- factor(database$type, levels = c("Peritumoral", "Tumor"))
        annot <- column_ha_type
    }

    ## Use gene names instead of ENSEMBL IDs
    mat <- SummarizedExperiment::assay(database)
    rownames(mat) <- rowData(database)$external_gene_name

    if (units == "log_TPM") mat <- log1p(mat)

    if (dim(mat)[1] > 100) fontsize <- 4
    if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) fontsize <- 5
    if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) fontsize <- 6
    if (dim(mat)[1] <= 20) fontsize <- 8

    h <- suppressMessages(Heatmap(mat[, , drop = FALSE],
                                  name = units,
                                  column_title = paste0("Expression in TCGA samples (", tumor, ")"),
                                  column_split = split_by,
                                  col = colorRamp2(seq(0, max(mat), length = 11),
                                                   legend_colors),
                                  clustering_method_rows = "ward.D",
                                  clustering_method_columns = "ward.D",
                                  cluster_rows = TRUE,
                                  show_column_names = FALSE,
                                  cluster_columns = TRUE,
                                  show_column_dend = FALSE,
                                  show_row_dend = FALSE,
                                  row_names_gp = gpar(fontsize = fontsize),
                                  heatmap_legend_param = legends_param,
                                  top_annotation = annot))

    if (return)
        return(mat)

    print(h)
}
