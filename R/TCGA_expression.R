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
#' If return = TRUE, the SummarizedExperiment with TPM expression data is
#' returned instead (Columns of RowData and colData are described in TCGA_TPM).
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' TCGA_expression(tumor = "LUAD",
#' genes = c("MAGEA1", "MAGEA3", "MAGEA4"),
#' units = "log_TPM")
TCGA_expression <- function(tumor = "all", genes = NULL,
                            units = "TPM", return = FALSE) {
  database <- TCGA_TPM

  if (!tumor %in% c(unique(sub(pattern = 'TCGA-', x = database$project_id, '')),
                    'all')) {
    stop("TCGA tumor code must be one of ",
         paste(c(unique(sub(pattern = 'TCGA-', x = database$project_id, ''))),
               ", "), "or all")
  }

  # Use only primary/metastatic tumors and normal peritumoral samples
  database <- database[, colData(database)$shortLetterCode %in% c("TP", "TM", "NT")]

  if (is.null(genes)) {
    database <- database[rowData(database)$external_gene_name %in% CT_genes$external_gene_name, ]
  }

  if (!is.null(genes)) {
    if (!all(genes %in% rowData(database)$external_gene_name)) {
      message("Check gene name(s)!\n")
      message(paste0(genes[!genes %in% rowData(database)$external_gene_name],
                     " is not in the database.\n"))
      genes <- genes[genes %in% rowData(database)$external_gene_name]
    }
    database <- database[rowData(database)$external_gene_name %in% genes, ]
  }

  if (tumor != "all") {
    database <- database[, colData(database)$project_id == paste0("TCGA-", tumor)]
  }

  if (tumor == "all") {
    database <- database[, colData(database)$shortLetterCode != "NT"]
  }

  colData(database)[colData(database)$definition == "Solid Tissue Normal",
                "definition"] <- "Peritumoral"
  colData(database)[colData(database)$definition == "Primary solid Tumor",
                "definition"] <- "Primary"
  df_col <- data.frame("barcode" = rownames(colData(database)),
                       "shortLetterCode" = colData(database)$shortLetterCode,
                       "tumor" = sub("TCGA-", x = colData(database)$project_id,
                                     replacement = ''),
                       "type" = factor(colData(database)$definition,
                                       levels = c("Peritumoral", "Primary",
                                                  "Metastatic")))

  df_col$tissue <- ifelse(df_col$shortLetterCode == "NT", "Peritumoral", "Tumor")
  rownames(df_col) <- rownames(colData(database))
  df_col <- df_col[order(df_col$tumor, df_col$tissue), ]

  legends_param <- list(
    labels_gp = gpar(col = "black", fontsize = 6),
    title_gp = gpar(col = "black", fontsize = 6),
    row_names_gp = gpar(fontsize = 4),
    annotation_name_side = "left")

  column_ha_type <- HeatmapAnnotation(
    Type = df_col$type,
    border = TRUE,
    col = list(Type = c("Metastatic" = "firebrick2", "Primary" = "salmon",
                        "Peritumoral" = "mediumseagreen")),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = legends_param)

  column_ha_tumor <- HeatmapAnnotation(
    Tumor = df_col$tumor,
    border = TRUE,
    col = list(Tumor = c("BRCA" = "pink", "COAD" = "black", "ESCA" = "green",
                         "HNSC" = "steelblue", "LUAD" = "cyan",
                         "LUSC" = "gray50", "SKCM" = "moccasin")),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = legends_param)

  split_order_by_cat <- factor(df_col$tissue,
                               levels = c("Peritumoral", "Tumor"))
  if (tumor == "all") {
    split_order_by_cat <- factor(df_col$tumor,
                                 levels = c("BRCA", "COAD", "ESCA", "HNSC",
                                            "LUAD", "LUSC", "SKCM"))
  }

  ## Use gene names instead of ENSEMBL IDs
  mat <- assay(database)
  rownames(mat) <- rowData(database)$external_gene_name
  name <- "TPM"
  title <- paste0("Gene expression in TCGA-", tumor)

  if (units == "log_TPM") {
    mat <- log1p(mat)
    name <- "log_TPM"
  }

  if (dim(mat)[1] > 100) { fontsize <- 4 }
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) { fontsize <- 5 }
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) { fontsize <- 6 }
  if (dim(mat)[1] <= 20) { fontsize <- 8 }

  annot <- column_ha_type
  if (tumor == "all") {
    annot <- column_ha_tumor
    title <- "Gene expression in TCGA samples"
  }


  h <- suppressMessages(Heatmap(mat[, rownames(df_col), drop = FALSE],
                                name = name,
                                column_title = title,
                                column_split = split_order_by_cat,
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
                                show_column_names = FALSE,
                                cluster_columns = TRUE,
                                show_column_dend = FALSE,
                                show_row_dend = FALSE,
                                row_names_gp = gpar(fontsize = fontsize),
                                heatmap_legend_param = legends_param,
                                top_annotation = annot))

  if (return == FALSE) {
    print(h)
  } else {
    mat
  }
}
