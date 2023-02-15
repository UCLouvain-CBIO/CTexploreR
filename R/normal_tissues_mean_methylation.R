#' Promoter methylation of Cancer-Testis genes in normal tissues
#'
#' @description Plots a heatmap of mean promoter methylation levels of
#'     Cancer-Testis (CT) genes in normal tissues. Methylation levels
#'     in tissues correspond to the mean methylation of CpGs located
#'     in range of 1000 pb upstream and 200 pb downstream from gene
#'     TSS.
#'
#' @param genes CT gene names.
#'
#' @param return `logical(1)`, `FALSE` by default. If `TRUE`, the
#'     function will return the methylation values in all samples
#'     instead of the heatmap.
#'
#' @param na.omit `logical(1)` specifying if genes with missing
#'     methylation values in some tissues should be removed (`TRUE` by
#'     default). Note that no gene clustering will be done when
#'     methylation values are missing.
#'
#' @return Heatmap of mean promoter methylation of Cancer-Testis (CT)
#'     genes in normal tissues. If return = TRUE, methylation values
#'     are returned instead.
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#'
#' @examples
#' normal_tissues_mean_methylation()
#' normal_tissues_mean_methylation(c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"))
#' normal_tissues_mean_methylation(c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"),
#'                                 na.omit = FALSE)
normal_tissues_mean_methylation <- function(genes = NULL, return = FALSE,
                                            na.omit = TRUE) {

    suppressMessages({
        database <- CTdata::CT_mean_methylation_in_tissues()
        CT_genes <- CTdata::CT_genes()
    })

    if (is.null(genes)) genes <- CT_genes$external_gene_name
    valid_gene_names <- unique(rownames(database))
    genes <- check_names(genes, valid_gene_names)
    database <- database[rownames(database) %in% genes]

    if (na.omit) {
        mat <- na.omit(assay(database))
        clustering_option <- TRUE
    } else {
        mat <- assay(database)
        clustering_option <- FALSE
    }

    if (dim(mat)[1] > 140 ) fontsize <- 3
    if (dim(mat)[1] > 100 & dim(mat)[1] <= 140) fontsize <- 4
    if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) fontsize <- 5
    if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) fontsize <- 6
    if (dim(mat)[1] <= 20) fontsize <- 8

    h <- Heatmap(mat,
                 column_title = 'Promoter mean methylation level by tissue',
                 name = 'Meth',
                 col = colorRamp2(c(1:100),
                                  colorRampPalette(c("moccasin","dodgerblue4"))
                                  (100)),
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

    if (return)
        return(mat)

    return(h)
}
