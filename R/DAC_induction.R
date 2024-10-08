#' Gene expression in cells treated or not by a demethylating agent
#'
#' @description Plots a heatmap of normalised gene counts
#'     (log-transformed) in a selection of cells treated or not by
#'     5-Aza-2'-Deoxycytidine (DAC), a demethylating agent.
#'
#' @param genes `character` naming the selected genes. The default
#'     value, `NULL`, takes all CT specific genes.
#'  
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#'
#' @param multimapping `logical(1)` defining whether to use
#'     multi-mapped gene expression dataset
#'     `CTdata::DAC_treated_cells_multimapping` or
#'     `DAC_treated_cells`. Default is `TRUE`.
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#'     gene normalised logcounts in all samples instead of the
#'     heatmap. Default is `FALSE`.
#'
#' @details
#'
#' RNAseq data from cells treated or not with 5-aza downloaded from
#' SRA.  (SRA references and information about cell lines and DAC
#' treatment are stored the colData of `DAC_treated_cells`).  Data was
#' processed using a standard RNAseq pipeline.
#' [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used
#' to align reads to grch38 genome.
#' [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html)
#' was used to assign reads to genes. Note that -M parameter was used
#' or not to allow or not counting multi-mapping reads.
#'
#' @return A heatmap of selected genes in cells treated or not by a
#'     demethylating agent. If `values_only` is `TRUE`, gene normalised
#'     logcounts are returned instead.
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
#'     multimapping = FALSE))
DAC_induction <- function(genes = NULL, multimapping = TRUE, 
                          include_CTP = FALSE,
                          values_only = FALSE) {
    
  suppressMessages({
    if (multimapping) {
      database <- CTdata::DAC_treated_cells_multimapping()
      } else {
        database <- CTdata::DAC_treated_cells()
        }
    })

    database <- subset_database(genes, database, include_CTP)

    mat <- assay(database, "log1p")
    rownames(mat) <- rowData(database)$external_gene_name

    df_col <- data.frame("cell" = colData(database)$cell,
                         "treatment" = colData(database)$treatment)
    rownames(df_col) <- rownames(colData(database))
    df_col <- df_col[order(df_col$cell, df_col$treatment), ]

    column_ha_cell <- HeatmapAnnotation(
        cell = df_col$cell,
        border = TRUE,
        col = list(cell = DAC_colors))

    column_ha_treatment <- HeatmapAnnotation(
        treatment = df_col$treatment,
        col = list(treatment = c("CTL" = "dodgerblue3", "DAC" = "firebrick1")),
        border = TRUE)

    fontsize <- set_fontsize(mat)

    h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
        name = "logCounts",
        column_title = "Gene expression in cells treated or not with 5-Aza",
        column_split = factor(df_col$cell),
        col = colorRamp2(seq(0, max(mat), length = 11), legend_colors),
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        clustering_method_rows = "ward.D",
        show_column_names = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = fontsize),
        top_annotation = c(column_ha_cell, column_ha_treatment))

    ifelse(values_only, return(mat), return(h))

}
