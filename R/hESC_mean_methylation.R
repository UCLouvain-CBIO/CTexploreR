#' Promoter methylation of any gene in hESC
#'
#' @description Plots a heatmap of mean promoter methylation levels of
#'     any genes in human embryonic cell lines. WGBS methylation data was 
#'     downloaded from Encode. Methylation levels in tissues correspond to 
#'     the mean methylation of CpGs located in range of 1000 pb upstream 
#'     and 200 pb downstream from gene TSS.
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
#' @return Heatmap of mean promoter methylation of any
#'     gene in hESC. If `values_only = TRUE`, a SummarizedExperiment cobtaining
#'     methylation values is returned instead.
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom SummarizedExperiment rowData<-
#' 
#'
#' @examples
#' \dontrun{
#' hESC_mean_methylation()
#' }
hESC_mean_methylation <- function(genes = NULL, 
                                  include_CTP = FALSE, 
                                  values_only = FALSE) {
  suppressMessages({
    database <- CTdata::mean_methylation_in_hESC()
  })
  
  rowData(database)$external_gene_name <- rownames(database)
  database <- subset_database(genes, database, include_CTP)
  
  mat <- assay(database)
  fontsize <- set_fontsize(mat)
  
  h <- Heatmap(mat,
               column_title = "Promoter mean methylation in human embryonic stem cells",
               name = "Meth",
               col = colorRamp2(seq_len(100),
                                colorRampPalette(c("moccasin", "dodgerblue4"))(100)),
               na_col = "gray80",
               cluster_rows = clustering_option,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_heatmap_legend = TRUE,
               show_row_dend = FALSE,
               row_names_gp = gpar(fontsize = fontsize),
               column_names_gp = gpar(fontsize = 12),
               column_names_rot = 0,
               column_names_side = "top",
               row_names_side = "right")
  
  ifelse(values_only, return(database), return(h))
}
