#' Promoter methylation of Cancer-Testis genes in normal tissues
#'
#' @description Plots a heatmap of mean promoter methylation levels of
#' Cancer-Testis (CT) genes in normal tissues. Methylation levels in
#' tissues correspond to the mean methylation of CpGs located in range of
#' 1000 pb upstream and 200 pb downstream from gene TSS.
#'
#' @param database CT_mean_methylation_in_tissues
#'
#' @param genes name of CT gene selected
#'
#' @param include_genes_with_missing_values Set to TRUE or FALSE to specify if
#' genes with missing methylation values in some tissues should be included
#' (set to FALSE by default). Note that no gene clustering will be done when
#' some methylation values are missing.
#'
#' @return heatmap of mean promoter methylation of Cancer-Testis (CT) genes
#' in normal tissues. Methylation values are returned invisibly.
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' normal_tissues_mean_methylation(CT_mean_methylation_in_tissues)
#' normal_tissues_mean_methylation(CT_mean_methylation_in_tissues,
#' c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"))
#' normal_tissues_mean_methylation(CT_mean_methylation_in_tissues,
#' c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"),
#' include_genes_with_missing_values = TRUE)
normal_tissues_mean_methylation <- function(database, genes = NULL,
                                            include_genes_with_missing_values = FALSE){

  if (missing(database)){
    stop("Database must be specified!")
  }

  if (!missing(database)){
    data <- database
  }

  if (!is.null(genes)) {
    if (!all(genes %in% rowData(data)$external_gene_name)) { ## Warning !
      message("Check gene name(s)!\n")
      message(paste0(genes[!genes %in% rowData(data)$external_gene_name], " is not in the database.\n"))
      genes <- genes[genes %in% rowData(data)$external_gene_name]
      stopifnot(length(genes) > 0)
    }
    data <- data[rowData(data)$external_gene_name %in% genes]
  }

  mat <- na.omit(assay(data))
  clustering_option <- TRUE

  if (include_genes_with_missing_values == TRUE){
    mat <- assay(data)
    clustering_option <- FALSE
  }
  if (dim(mat)[1] > 140 ){ fontsize <- 3 }
  if (dim(mat)[1] > 100 & dim(mat)[1] <= 140){ fontsize <- 4 }
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100){ fontsize <- 5 }
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50){ fontsize <- 6 }
  if (dim(mat)[1] <= 20) { fontsize <- 8 }

  h <- Heatmap(mat,
               column_title = 'Promoter mean methylation level by tissue',
               name = 'Methylation',
               col = colorRamp2(c(1:100),
                                colorRampPalette(c("moccasin","dodgerblue4"))(100)),
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

  print(h)
  invisible(mat)
}

