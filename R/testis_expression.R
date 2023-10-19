#' Gene expression in testis cells
#'
#' @description
#'
#' Plots a heatmap of genes expression in the different types of testis cells,
#' using scRNAseq data from "The adult human testis transcriptional cell atlas"
#' (Guo et al. 2018)
#'
#' @param cells `character` defining the testis cell types to be plotted.
#' Can be "germ_cells", "somatic_cells", "all" (default), or any or a
#' combination of "SSC", "Spermatogonia", "Early_spermatocyte",
#' "Late_spermatocyte", "Round_spermatid", "Elongated_spermatid", "Sperm1",
#' "Sperm2", "Macrophage", "Endothelial", "Myoid", "Sertoli", "Leydig".
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT genes.
#'
#' @param scale_lims `vector of length 2` setting the lower and upper limits
#' of the heatmap colorbar. By default, the lower limit is 0, and the upper
#' limit corresponds to the third quartile of the logcounts values.
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#'     SingleCellExperiment instead of the heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected CT genes expression in single cells from adult
#'     testis. If `values_only = TRUE`, a SingleCellExperiment instead of the
#'     heatmap is returned instead.
#'
#' @export
#'
#' @importFrom SingleCellExperiment logcounts colData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom MatrixGenerics rowMaxs
#'
#' @examples
#' \dontrun{
#' testis_expression(cells = "germ_cells",
#'                   genes = c("MAGEA1", "MAGEA3", "MAGEA4"))
#' }
testis_expression <- function(
    cells = c("all", "germ_cells", "somatic_cells", "SSC", "Spermatogonia",
              "Early_spermatocyte", "Late_spermatocyte", "Round_spermatid", 
              "Elongated_spermatid", "Sperm1", "Sperm2", "Macrophage", 
              "Endothelial", "Myoid", "Sertoli", "Leydig"), 
    genes = NULL, scale_lims = NULL, values_only = FALSE) {
  
    suppressMessages({
        database <- CTdata::testis_sce()
        CT_genes <- CTdata::CT_genes()
    })
  
    cells <- match.arg(cells, several.ok = TRUE)
  
    germ_cells <- c("SSC", "Spermatogonia", "Early_spermatocyte",
                  "Late_spermatocyte", "Round_spermatid", "Elongated_spermatid",
                  "Sperm1", "Sperm2")
  
    somatic_cells <- c("Macrophage", "Endothelial",
                     "Myoid", "Sertoli", "Leydig")
  
    if ("all" %in% cells) { 
      cells <- c(germ_cells, somatic_cells) 
    } else if (length(cells) > 1) {
      stopifnot("`cells` parameter can be set to 'all', 'germ_cells', 
        'somatic_cells', or any combination of testis cell type" 
                = all(cells %in% c(germ_cells, somatic_cells)))
    } else if (!cells %in% c(germ_cells, somatic_cells)) {
      cells <- switch(cells, 
                      "germ_cells" = germ_cells, 
                      "somatic_cells" = somatic_cells)
    }
    
    database <- subset_database(genes, database[, database$type %in% cells])

    mat <- SingleCellExperiment::logcounts(database)

    df_col <- data.frame(clusters = colData(database)$clusters,
                         type = colData(database)$type,
                         Donor = colData(database)$Donor)
    rownames(df_col) <- colnames(database)
    df_col <- df_col[order(df_col$type), ]
    df_col$lineage <- "Germ cells"
    df_col$lineage[df_col$type %in% somatic_cells] <- "Somatic cells"

    column_ha_type <- HeatmapAnnotation(
        type = df_col$type,
        border = TRUE,
        col = list(type = testis_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    column_ha_lineage <- HeatmapAnnotation(
        lineage = df_col$lineage,
        border = TRUE,
        col = list(lineage = c(
            "Germ cells" = "salmon",
            "Somatic cells" = "cyan4")),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    fontsize <- set_fontsize(mat)

    if (is.null(scale_lims)) scale_lims <- c(0, quantile(rowMaxs(mat), 0.75))

    if (any(database$type %in% germ_cells) &
        any(database$type %in% somatic_cells)) {
        top_annot <- c(column_ha_lineage, column_ha_type)
    } else {
        top_annot <- column_ha_type
    }

    h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
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

    if (values_only) {
        return(database)
    }

    return(h)
}


