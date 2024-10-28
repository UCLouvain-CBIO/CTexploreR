#' Promoter methylation of any gene in fetal germ cells
#'
#' @description Plots a heatmap of mean promoter methylation levels of
#'     any genes in fetal germ cells, using WGSB data from "Dissecting the 
#'     epigenomic dynamics of human fetal germ cell development at single-cell 
#'     resolution" (Li et al. 2021). Methylation levels in tissues correspond 
#'     to the mean methylation of CpGs located in range of 1000 pb upstream and 
#'     500 pb downstream from gene TSS.
#'
#' @param genes `character` naming the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'     
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#'
#' @param values_only `logical(1)`, `FALSE` by default. If `TRUE`, the
#'     function will return the methylation values in all samples
#'     instead of the heatmap.
#'
#' @return Heatmap of mean promoter methylation of any gene in normal tissues. 
#' If `values_only = TRUE`, a SummarizeExperiment with methylation values is
#'  returned instead.
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar unit
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment rowData<- assay
#' 
#'
#' @examples
#' fetal_germcells_mean_methylation()
#' fetal_germcells_mean_methylation(c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEC2"))
fetal_germcells_mean_methylation <- function(genes = NULL, 
                                            include_CTP = FALSE, 
                                            values_only = FALSE) {
    suppressMessages({
        database <- CTdata::mean_methylation_in_FGC()
    })

    database <- subset_database(genes, database, include_CTP)
    
    df_col <- colData(database)
    
    df_col <- df_col[order(df_col$sex_time),] 
    
    column_ha_time = HeatmapAnnotation(
      time_week = df_col$time_week,
      col = list(time_week = fetal_time_colors),   
      border = FALSE, 
      simple_anno_size = unit(0.3,"cm"), 
      annotation_legend_param = legends_param,
      annotation_name_gp = gpar(fontsize = 6))
    
    column_ha_sex = HeatmapAnnotation(
      sex = df_col$sex,
      col = list(sex = sex_colors),   
      border = FALSE, 
      simple_anno_size = unit(0.3,"cm"), 
      annotation_legend_param = legends_param,
      annotation_name_gp = gpar(fontsize = 6))
    
    column_ha_type = HeatmapAnnotation(
      type = df_col$type,
      col = list(type = fetal_type_colors),
      border = FALSE, 
      simple_anno_size = unit(0.3,"cm"), 
      annotation_legend_param = legends_param,
      annotation_name_gp = gpar(fontsize = 6))

    mat <- assay(database[, df_col$sample])
    fontsize <- set_fontsize(mat)

    h <- Heatmap(mat,
        column_title = "Promoter mean methylation level by tissue",
        name = "Meth",
        col = colorRamp2(seq_len(100),
                         colorRampPalette(c("moccasin", "dodgerblue4"))(100)),
        na_col = "gray80",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_heatmap_legend = TRUE,
        show_row_dend = FALSE,
        column_split = df_col$sex_type,
        column_gap = unit(c(1), "mm"),
        top_annotation = c(column_ha_type,
          column_ha_sex,
          column_ha_time,
          gap = unit(0.1, "mm")), 
        row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = 0),
        column_names_side = "bottom",
        row_names_side = "right")

    ifelse(values_only, return(database[, df_col$sample]), return(h))
}
