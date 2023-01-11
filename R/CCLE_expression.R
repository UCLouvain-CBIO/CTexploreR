#' Gene expression in CCLE Tumors
#'
#' @description Plots an expression heatmap of genes in CCLE tumor cell lines.
#'
#' @param type `character` describing the tumor cell line(s) type to
#' be plotted. Allowed cell lines are "Ovarian", "Leukemia", "Colorectal",
#' Skin", "Lung", "Bladder", "Kidney", "Breast", "Pancreatic", "Myeloma",
#' "Brain", "Sarcoma", "Lymphoma", "Bone", "Neuroblastoma", "Gastric",
#' "Uterine", "Head_and_Neck", "Bile_Duct" and "Esophageal".
#'
#' @param genes Genes selected (all CT genes by default)
#'
#' @param units Expression values unit.
#' Can be "TPM" (default) or "log_TPM" (log(TPM + 1))
#'
#' @param database CCLE_data
#'
#' @return A heatmap of selected genes in CCLE cell lines from specified type.
#' Expression values are invisibly returned.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' CCLE_expression(genes = c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEA6", "MAGEA10"),
#'                 type = c("Skin", "Lung"), units = "log_TPM")
#' CCLE_expression(units = "log_TPM")
CCLE_expression <- function(genes = NULL, units = "TPM", type = NULL,
                            database = CCLE_data) {

  CCLE <- database

  if (!is.null(type)) {
    if (!all(type %in% colData(CCLE)$type)) {
      message("Tumor type(s) must be one of: ",
              paste0(unique(colData(CCLE)$type), ' '))
      return(invisible(NA))
    }
    CCLE <- CCLE[, colData(CCLE)$type %in% type]
  }

  if (is.null(genes)) {
    CCLE <- CCLE[rowData(CCLE_data)$external_gene_name %in% CT_genes$external_gene_name, ]
  }

  if (!is.null(genes)) {
    if (!all(genes %in% rowData(CCLE)$external_gene_name)) {
      message("Check gene name(s)!\n")
      message(paste0(genes[!genes %in% rowData(CCLE)$external_gene_name],
                     " is not in the database.\n"))
    }
    CCLE <- CCLE[rowData(CCLE)$external_gene_name %in% genes, ]
  }

  mat <- assay(CCLE)
  rownames(mat) <- rowData(CCLE)$external_gene_name

  name <- "TPM"
  if (units == "log_TPM") {
    mat <- log1p(mat)
    name <- "log_TPM"
  }

  legends_param <- list(
    labels_gp = gpar(col = "black", fontsize = 6),
    title_gp = gpar(col = "black", fontsize = 6),
    row_names_gp = gpar(fontsize = 4),
    annotation_name_side = "left")

  df_col <- data.frame("cell_line" = colData(CCLE)$cell_line_name,
                       "type" = as.factor(colData(CCLE)$type))
  rownames(df_col) <- rownames(colData(CCLE))
  df_col <- df_col[order(df_col$type), ]

  set.seed(1)

  column_ha_type <- HeatmapAnnotation(
    type = df_col$type,
    border = TRUE,
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = legends_param)

  if (dim(mat)[1] > 100) { fontsize <- 4 }
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) { fontsize <- 5 }
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) { fontsize <- 6 }
  if (dim(mat)[1] <= 20) { fontsize <- 8 }

  if (length(type) <= 5) { label_fontsize <- 6 }
  if (length(type) > 5 & length(type) < 10) { label_fontsize <- 4 }
  if (length(type) >= 10 | is.null(type)) { label_fontsize <- 0 }

  h <- suppressMessages(Heatmap(mat[, rownames(df_col), drop = FALSE],
                                name = name,
                                column_split = factor(df_col$type),
                                col = colorRamp2(seq(0, max(mat), length = 11),
                                                 c("#5E4FA2", "#3288BD",
                                                   "#66C2A5", "#ABDDA4",
                                                   "#E6F598", "#FFFFBF",
                                                   "#FEE08B", "#FDAE61",
                                                   "#F46D43", "#D53E4F",
                                                   "#9E0142")),
                                clustering_method_rows = "ward.D",
                                clustering_method_columns = "ward.D",
                                cluster_rows = TRUE,
                                show_row_dend = FALSE,
                                show_column_names = FALSE,
                                cluster_columns = TRUE,
                                show_column_dend = FALSE,
                                column_title_gp = gpar(fontsize = label_fontsize),
                                row_names_gp = gpar(fontsize = fontsize),
                                heatmap_legend_param = legends_param,
                                top_annotation = c(column_ha_type)))

  print(h)
  invisible(mat)
}
