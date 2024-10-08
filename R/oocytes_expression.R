#' Gene expression in oocytes
#'
#' @description
#'
#' Plots a heatmap of genes expression in oocytes,
#' using scRNAseq data from "Decoding dynamic epigenetic landscapes in human 
#' oocytes using single-cell multi-omics sequencing"
#' (Yan et al. Cell Stem Cell 2021)
#'
#' @param genes `character` naming the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'     
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#' 
#' @param ncells_max `integer(1)` Sets the number of each cell type to 
#' represent on the heatmap (these cells will be randomly selected among each 
#' cell type) (set to 200 by default). If `NULL`, all cells are displayed.
#'
#' @param scale_lims `vector of length 2` setting the lower and upper limits
#' of the heatmap colorbar. By default, the lower limit is 0, and the upper
#' limit corresponds to the third quartile of the logcounts values.
#'
#' @param values_only `logical(1)`. If `TRUE`, the function will return the
#'     SingleCellExperiment instead of the heatmap. Default is `FALSE`.
#'
#' @return A heatmap of selected CT genes expression in single cells from human
#' oocytes. If `values_only = TRUE`, a SingleCellExperiment is returned instead.
#'
#' @export
#'
#' @importFrom SingleCellExperiment logcounts colData rowData
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom MatrixGenerics rowMaxs
#'
#' @examples
#' \dontrun{
#' oocytes_expression(include_CTP = FALSE, ncells_max = 100,
#'                   values_only = FALSE)
#' }
oocytes_expression <- function(
    genes = NULL, include_CTP = FALSE, 
    ncells_max = 200,
    scale_lims = NULL, values_only = FALSE) {
  
    suppressMessages({
        database <- CTdata::oocytes_sce()
        rowData(database)$external_gene_name <- rownames(database)
    })
    
    database <- subset_database(genes, database, include_CTP)

    database <- database[, database$type %in% c("Growing oocytes", 
                                                "Fully grown oocytes", 
                                                "Metaphase I", 
                                                "Metaphase II" )]
    database$type <- droplevels(database$type)

    #Randomly pick ncells_max cells of each cell type
    if (!is.null(ncells_max)) {
      set.seed(123)
      selected_cells <- sample(
        colnames(database[, database$type == levels(database$type)[1]]), 
        size = min(ncells_max, length(colnames(
          database[, database$type == levels(database$type)[1]]))), 
        replace = FALSE)
        
        for (type in levels(database$type)[-1]){
          tmp <- sample(
            colnames(database[, database$type == type]), 
            size = min(ncells_max, length(colnames(
              database[, database$type == type]))), 
            replace = FALSE)
          selected_cells <- c(selected_cells, tmp)
        }
      database <- database[, selected_cells]
    }
    
    mat <- SingleCellExperiment::logcounts(database)

    df_col <- data.frame(sex = database$sex,
                         type = database$type)
    rownames(df_col) <- colnames(database)
    df_col <- df_col[order(df_col$type), ]

    column_ha_type <- HeatmapAnnotation(
        type = df_col$type,
        border = TRUE,
        col = list(type = oocytes_colors),
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = legends_param)

    fontsize <- set_fontsize(mat)

    if (is.null(scale_lims)) scale_lims <- c(0, quantile(rowMaxs(mat), 0.75))

    h <- Heatmap(mat[, rownames(df_col), drop = FALSE],
        name = "logCounts",
        column_title = "Expression in oocytes (scRNAseq)",
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
        top_annotation = column_ha_type,
        heatmap_legend_param = legends_param)

    ifelse(values_only, return(database), return(h))
}

