#' Gene expression in testis cells
#'
#' @description
#'
#' Plots a heatmap of genes expression in the different types of testis cells, 
#' using scRNAseq data from "The adult human testis transcriptional cell atlas" 
#' (Guo et al. 2018)
#'
#' @param cells `character` defining the testis cell types to be plotted. 
#' Can be "germ_cells", "somatic_cells", "all" (default), or any or a combination 
#' of "SSC", "Spermatogonia", "Early_spermatocyte", "Late_spermatocyte", 
#' "Round_spermatid", "Elongated_spermatid", "Sperm1", "Sperm2", "Macrophage", 
#' "Endothelial", "Myoid", "Sertoli", "Leydig".
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT genes.
#'     
#' @param scale_lims `vector of length 2` setting the lower and upper limits of 
#'    the heatmap colorbar. By default, the lower limit is 0, and the upper limit
#'    corresponds to the third quartile of the logcounts values.
#'
#' @param return `logical(1)`. If `TRUE`, the function will return the
#'     SingleCellExperiment instead of the heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected CT genes expression in single cells from adult 
#'     testis. If `return = TRUE`, a SingleCellExperiment instead of the heatmap 
#'     is returned instead.
#'
#' @export
#'
#' @importFrom SingleCellExperiment logcounts colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom Biobase rowMax
#'
#' @examples
#'
#' testis_expression(cells = "germ_cells", genes = c("MAGEA1", "MAGEA3", "MAGEA4"))
#' testis_expression(cells = "all", scale_lims = c(0, 4))
testis_expression <- function(cells = "all", genes = NULL, 
                              scale_lims = NULL, return = FALSE) {
    suppressMessages({
      database <- CTdata::testis_sce()
      CT_genes <- CTdata::CT_genes()
    })
    
    germ_cells <- c("SSC", "Spermatogonia", "Early_spermatocyte", 
                  "Late_spermatocyte","Round_spermatid", "Elongated_spermatid",
                  "Sperm1", "Sperm2")
  
    somatic_cells <- c("Macrophage", "Endothelial", "Myoid", "Sertoli", "Leydig")
  
    valid_cell_names <- c("germ_cells", "somatic_cells", "all", germ_cells, 
                          somatic_cells)
    cells <- check_names(cells, valid_cell_names)
    
    if (all(cells %in% c("germ_cells"))) {
      database <- database[, database$type %in% germ_cells]
    } else if ((all(cells %in% c("somatic_cells")))) {
      database <- database[, database$type %in% somatic_cells]
    } else if (all(cells %in% c(germ_cells, somatic_cells))) {
      database <- database[, database$type %in% cells]
    }

    if (is.null(genes)) genes <- CT_genes$external_gene_name
    valid_gene_names <- unique(rownames(database))
    genes <- check_names(genes, valid_gene_names)
    database <- database[genes, ]
  
    mat <- SingleCellExperiment::logcounts(database)

    legends_param <- list(
        labels_gp = gpar(col = "black", fontsize = 6),
        title_gp = gpar(col = "black", fontsize = 6),
        row_names_gp = gpar(fontsize = 4),
        annotation_name_side = "left")
    
    df_col <- data.frame(clusters = colData(database)$clusters,
                         type = colData(database)$type,
                         Donor = colData(database)$Donor)
    rownames(df_col) <- colnames(database)
    df_col <- df_col[order(df_col$type),]
    df_col$lineage <- "Germ cells"
    df_col$lineage[df_col$type %in% somatic_cells] <- "Somatic cells"
    
    column_ha_type = HeatmapAnnotation(
      type = df_col$type,
      border = TRUE,
      col = list(type = testis_colors),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_legend_param = legends_param)
    
    column_ha_lineage = HeatmapAnnotation(
      lineage = df_col$lineage,
      border = TRUE,
      col = list(lineage = c("Germ cells" = "salmon", "Somatic cells" = "cyan4")),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_legend_param = legends_param)
      
    if (dim(mat)[1] > 100) fontsize <- 4
    if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) fontsize <- 5
    if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) fontsize <- 6
    if (dim(mat)[1] <= 20) fontsize <- 8
    
    if (is.null(scale_lims)) scale_lims <- c(0, quantile(rowMax(mat), 0.75))
    
    if (any(database$type %in% germ_cells) & 
        any(database$type %in% somatic_cells)) {
      top_annot <- c(column_ha_lineage, column_ha_type)
    } else { top_annot <- column_ha_type
    }
    
    h <- Heatmap(mat[genes, rownames(df_col), drop = FALSE],
            name = "logCounts",
            column_title = "Expression in testis cells (scRNAseq)",
            column_split = df_col$type,
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
            top_annotation = top_annot,
            heatmap_legend_param = legends_param)

    if (return)
        return(database)

    return(h)
}

