#' Gene expression in cells treated or not by a demethylating agent
#'
#' @description Plots a heatmap of normalised gene counts (log-transformed)
#' in a selection of cells treated or not by 5-Aza-2′-Deoxycytidine (DAC),
#' a demethylating agent.
#'
#' @param database Can be `DAC_treated_cells` or` DAC_treated_cells_multimapping`,
#' depending if returned expression values should take into account or not
#' multi-mapped reads.
#'
#' @param genes Genes selected (all CT genes by default)
#'
#' @details
#' RNAseq data from cells treated or not with 5-aza downloaded from SRA.
#' Data was processed using a standard RNAseq pipeline.
#' [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
#' reads to grch38 genome.
#' [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was used
#' to assign reads to genes. Note that -M parameter was used or not to allow or not
#' counting multi-mapping reads.
#'
#' @return A heatmap of selected genes in cells treated or not by a demethylating
#' agent. Gene normalised logcounts are invisibly returned.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' DAC_induction(database = DAC_treated_cells, genes = c("MAGEA1", "MAGEA3",
#' "MAGEA4", "MAGEA6", "CTAG1A"))
#' DAC_induction(database = DAC_treated_cells_multimapping,
#' genes = c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEA6", "CTAG1A"))
DAC_induction <- function(database, genes = NULL) {

  if (missing(database)) {
    stop("Database must be specified!")
  }

  if (!missing(database)) {
    data <- database
  }

  if (!is.null(genes)) {
    if (!all(genes %in% rowData(data)$external_gene_name)) {
      message("Check gene name(s)!\n")
      message(paste0(genes[!genes %in% rowData(data)$external_gene_name],
                     " is not in the database.\n"))
      genes <- genes[genes %in% rowData(data)$external_gene_name]
      stopifnot(!is_empty(genes))
    }
    data <- data[rowData(data)$external_gene_name %in% genes]
  }

  if (is.null(genes)) {
    data <- data[rowData(data)$external_gene_name %in% CT_genes$external_gene_name, ]
  }

  mat <- assay(data, "log1p")
  rownames(mat) <- rowData(data)$external_gene_name

  set.seed(1)
  df_col <- data.frame("cell" = colData(data)$cell,
                       "treatment" = colData(data)$treatment)
  rownames(df_col) <- rownames(colData(data))
  df_col <- df_col[order(df_col$cell, df_col$treatment), ]

  column_ha_cell <- HeatmapAnnotation(cell = df_col$cell,border = TRUE)

  column_ha_treatment <- HeatmapAnnotation(
    treatment = df_col$treatment,
    col = list(treatment = c("CTL" = "cyan", "DAC" = "firebrick1")),
    border = TRUE)

  if (dim(mat)[1] > 100) { fontsize <- 4 }
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) { fontsize <- 5 }
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) { fontsize <- 6 }
  if (dim(mat)[1] <= 20) { fontsize <- 8 }

  h <- suppressMessages(Heatmap(mat[, rownames(df_col), drop = FALSE],
                                name = "logCounts",
                                column_split = factor(df_col$cell),
                                col = colorRamp2(seq(0, max(mat), length = 11),
                                                 c("#5E4FA2", "#3288BD",
                                                   "#66C2A5", "#ABDDA4",
                                                   "#E6F598", "#FFFFBF",
                                                   "#FEE08B", "#FDAE61",
                                                   "#F46D43", "#D53E4F",
                                                   "#9E0142")),

                                cluster_rows = TRUE,
                                show_row_dend = FALSE,
                                clustering_method_rows = "ward.D",
                                show_column_names = FALSE,
                                cluster_columns = FALSE,
                                row_names_gp = gpar(fontsize = fontsize),
                                column_title_gp = gpar(fontsize = 0),
                                top_annotation = c(column_ha_cell,
                                                   column_ha_treatment)))
  print(h)
  invisible(mat)
}

