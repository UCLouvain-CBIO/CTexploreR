#' Gene expression in cells treated or not by a demethylating agent
#'
#' @description Plots a heatmap of normalised gene counts (log-transformed)
#' in a selection of cells treated or not by 5-Aza-2′-Deoxycytidine (DAC),
#' a demethylating agent.
#'
#' @param genes Genes selected (all CT genes by default)
#'
#' @param return Boolean (FALSE by default). If set to TRUE, the function will
#' return the gene normalised logcounts in all samples instead of the heatmap.
#'
#' @details
#' RNAseq data from cells treated or not with 5-aza downloaded from SRA.
#' (SRA references and information about cell lines and DAC treatment are stored
#' the colData of `DAC_treated_cells`).
#' Data was processed using a standard RNAseq pipeline.
#' [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
#' reads to grch38 genome.
#' [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was used
#' to assign reads to genes. Note that -M parameter was used or not to allow or not
#' counting multi-mapping reads.
#'
#' @return A heatmap of selected genes in cells treated or not by a demethylating
#' agent. If return = TRUE, gene normalised logcounts are returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData colData assay
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' DAC_induction(genes = c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEA6", "CTAG1A"))
#' DAC_induction(genes = c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEA6", "CTAG1A",
#' multimapping = FALSE))
DAC_induction <- function(genes = NULL, multimapping = TRUE, return = FALSE) {

  if (multimapping == TRUE) {
    database <- DAC_treated_cells_multimapping
  } else {
    database <- DAC_treated_cells
  }

  if (is.null(genes)) genes <- CT_genes$external_gene_name
  valid_gene_names <- unique(rowData(database)$external_gene_name)
  genes <- check_names(genes, valid_gene_names)
  database <- database[rowData(database)$external_gene_name %in% genes, ]

  mat <- assay(database, "log1p")
  rownames(mat) <- rowData(database)$external_gene_name

  set.seed(1)
  df_col <- data.frame("cell" = colData(database)$cell,
                       "treatment" = colData(database)$treatment)
  rownames(df_col) <- rownames(colData(database))
  df_col <- df_col[order(df_col$cell, df_col$treatment), ]

  column_ha_cell <- HeatmapAnnotation(cell = df_col$cell,
                                      border = TRUE)

  column_ha_treatment <- HeatmapAnnotation(
    treatment = df_col$treatment,
    col = list(treatment = c("CTL" = "cyan", "DAC" = "firebrick1")),
    border = TRUE)

  if (dim(mat)[1] > 100) fontsize <- 4
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) fontsize <- 5
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) fontsize <- 6
  if (dim(mat)[1] <= 20) fontsize <- 8

  h <- suppressMessages(Heatmap(mat[, rownames(df_col), drop = FALSE],
                                name = "logCounts",
                                column_title = "Gene expression in cells treated or not with 5-Aza",
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
                                top_annotation = c(column_ha_cell,
                                                   column_ha_treatment)))


  if (return == FALSE) {
    print(h)
  } else {
    mat
  }

}

