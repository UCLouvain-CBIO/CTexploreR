#' Gene expression in human embryos
#'
#' @description
#'
#' Plots a heatmap of genes expression in human early embryos, 
#' from "Petropoulos" scRNAseq dataset ("Single-Cell RNA-Seq Reveals Lineage 
#' and X Chromosome Dynamics in Human Preimplantation Embryos". 
#' Petropoulos et al., Cell 2016) or from "Zhu" scRNAseq dataset ("Single-cell 
#' DNA methylome sequencing of human preimplantation embryos". Zhu et al. 
#' Nat genetics 2018)
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'     
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#'
#' @param scale_lims `vector of length 2` setting the lower and upper limits
#' of the heatmap colorbar. By default, the lower limit is 0, and the upper
#' limit corresponds to the third quartile of the logcounts values.
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#'     SingleCellExperiment instead of the heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected CT genes expression in single cells from 
#' embryos. If `values_only = TRUE`, a SingleCellExperiment is returned instead.
#'
#' @export
#'
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom MatrixGenerics rowMaxs
#' 
#' @examples
#' \dontrun{
#' embryo_expression(dataset = "Petropoulos", include_CTP = FALSE)
#' embryo_expression(dataset = "Zhu", include_CTP = FALSE)
#' }
embryo_expression <- function(
    dataset = c("Petropoulos", "Zhu"),
    genes = NULL, include_CTP = FALSE, 
    scale_lims = NULL, values_only = FALSE) {
  
  if (dataset == "Petropoulos") {
    suppressMessages({
      database <- CTdata::embryo_sce_Petropoulos()
      rowData(database)$external_gene_name <- rownames(database)
      
      df_col <- data.frame(stage = database$day_stage,
                           sex = database$sex,
                           genotype = database$genotype,
                           day = database$day)
      rownames(df_col) <- colnames(database)
      df_col <- df_col[order(df_col$stage, df_col$sex), ]
      
      column_ha_stage <- HeatmapAnnotation(
        stage = df_col$stage,
        border = TRUE,
        col = list(stage = Petropoulos_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)
      
      column_ha_sex <- HeatmapAnnotation(
        sex = df_col$sex,
        border = TRUE,
        col = list(sex = sex_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)
      
      column_ha_genotype <- HeatmapAnnotation(
        genotype = df_col$genotype,
        border = TRUE,
        col = list(genotype = genotype_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)
      
      units <- "log_RPKM"
    })
  }
  
  if (dataset == "Zhu") {
    suppressMessages({
      database <- CTdata::embryo_sce_Zhu()
      rowData(database)$external_gene_name <- rownames(database)
      df_col <- data.frame(stage = database$stage,
                           sex = database$sex,
                           genotype = database$genotype)
      rownames(df_col) <- colnames(database)
      df_col <- df_col[order(df_col$stage, df_col$sex), ]
      
      column_ha_stage <- HeatmapAnnotation(
        stage = df_col$stage,
        border = TRUE,
        col = list(stage = c("blastocyst" = "gray")),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)
      
      column_ha_sex <- HeatmapAnnotation(
        sex = df_col$sex,
        border = TRUE,
        col = list(sex = sex_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)
      
      column_ha_genotype <- HeatmapAnnotation(
        genotype = df_col$genotype,
        border = TRUE,
        col = list(genotype = genotype_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)
      
      units <- "log_FPKM"
    })
  }
  
  database <- subset_database(genes, database, include_CTP)
  mat <- log1p(assay(database))
  
  fontsize <- set_fontsize(mat)
  if (is.null(scale_lims)) scale_lims <- c(0, quantile(rowMaxs(mat), 0.75))
  
  h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
               name = units,
               column_title = paste0("Expression in early embryos (scRNAseq, ", 
                                     dataset, " dataset)"),
               column_split = df_col$stage,
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
               top_annotation = c(column_ha_sex, column_ha_stage),
               heatmap_legend_param = legends_param)
  
  ifelse(values_only, return(database), return(h))
}





