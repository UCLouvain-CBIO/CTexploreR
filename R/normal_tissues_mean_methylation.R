#' Promoter methylation of any gene in normal tissues
#'
#' @description Plots a heatmap of mean promoter methylation levels of
#'     any genes in normal tissues. Methylation levels
#'     in tissues correspond to the mean methylation of CpGs located
#'     in range of 1000 pb upstream and 200 pb downstream from gene
#'     TSS.
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
#' @param na.omit `logical(1)` specifying if genes with missing
#'     methylation values in some tissues should be removed (`TRUE` by
#'     default). Note that no gene clustering will be done when
#'     methylation values are missing.
#'
#' @return Heatmap of mean promoter methylation of any
#'     gene in normal tissues. If `values_only = TRUE`, methylation values
#'     are returned instead.
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
#' normal_tissues_mean_methylation()
#' normal_tissues_mean_methylation(c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"))
#' normal_tissues_mean_methylation(c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"),
#'     na.omit = FALSE)
normal_tissues_mean_methylation <- function(genes = NULL, 
                                            include_CTP = FALSE, 
                                            values_only = FALSE,
                                            na.omit = TRUE) {
    suppressMessages({
        database <- CTdata::mean_methylation_in_tissues()
    })

    rowData(database)$external_gene_name <- rownames(database)
    database <- subset_database(genes, database, include_CTP)
    
    if (na.omit) {
        mat <- na.omit(assay(database))
        clustering_option <- TRUE
    } else {
        mat <- assay(database)
        clustering_option <- FALSE
    }
    
    fontsize <- set_fontsize(mat)

    h <- Heatmap(mat,
        column_title = "Promoter mean methylation level by tissue",
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
        column_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        row_names_side = "right")

    ifelse(values_only, return(mat), return(h))
}
