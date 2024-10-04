#' Promoter methylation of any gene in early embryos
#'
#' @description Plots a heatmap of mean promoter methylation levels of
#'     any genes in early embryos, using WGSB data from ("Single-cell 
#'     DNA methylome sequencing of human preimplantation embryos". Zhu et al. 
#'     Nat genetics 2018). Methylation levels in tissues correspond 
#'     to the mean methylation of CpGs located in range of 1000 pb upstream and 
#'     500 pb downstream from gene TSS.
#'
#' @param genes `character` naming the selected genes. The default
#'     value, `NULL`, takes all CT (specific) genes.
#'     
#' @param cells `character` defining the cell types to be plotted.
#' Can be "GV Oocyte", "MII Oocyte", "Sperm", "Zygote", "2-cell", 
#' "4-cell", "8-cell", "Morula", "Blastocyst", "Post-implantation".
#'     
#' @param include_CTP `logical(1)` If `TRUE`, CTP genes are included.
#' (`FALSE` by default).
#'
#' @param values_only `logical(1)`, `FALSE` by default. If `TRUE`, the
#'     function will return the methylation values in all samples
#'     instead of the heatmap.
#'
#' @return Heatmap of mean promoter methylation of any gene in embryos. 
#' If `values_only = TRUE`, a RangedSummarizedExperiment with methylation values 
#' is returned instead.
#'
#' @export
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom SummarizedExperiment rowData<- assay
#' 
#'
#' @examples
#' embryos_mean_methylation()
#' embryos_mean_methylation(c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEC2", "MAGEB16),
#' stage = c( "MII Oocyte", "Sperm", "Zygote", "2-cell", "4-cell", "8-cell", 
#' "Morula"))
embryos_mean_methylation <- function(genes = NULL, 
                                     stage = c("GV Oocyte", "MII Oocyte", 
                                               "Sperm", "Zygote", "2-cell",
                                               "4-cell", "8-cell", "Morula", 
                                               "Blastocyst", 
                                               "Post-implantation"),
                                     include_CTP = FALSE, 
                                     values_only = FALSE) {
    suppressMessages({
        #database <- CTdata::mean_methylation_in_FGC()
        load("/storage/research/dduv/cbio-lg/cluster/Packages/CTdata/eh_data/mean_methylation_in_embryo.rda")
        database <- mean_methylation_in_embryo
    })

    database <- subset_database(genes, database, include_CTP)
    stage <- match.arg(stage, several.ok = TRUE)
    database <- database[, database$Stage %in% stage]
    df_col <- colData(database)
    df_col <- df_col[order(df_col$Stage, df_col$`Genotype of the embryo`),] 
    
    column_ha_stage = HeatmapAnnotation(
      stage = df_col$Stage,
      col = list(stage = embryo_stage_col),   
      border = FALSE, 
      simple_anno_size = unit(0.3,"cm"), 
      annotation_legend_param = legends_param,
      annotation_name_gp = gpar(fontsize = 6))
    
    column_ha_genotype = HeatmapAnnotation(
      genotype = df_col$`Genotype of the embryo`,
      col = list(genotype = genotype_colors),   
      border = FALSE, 
      simple_anno_size = unit(0.3,"cm"), 
      annotation_legend_param = legends_param,
      annotation_name_gp = gpar(fontsize = 6))

    mat <- SummarizedExperiment::assay(database[, df_col$Sample_Name])
    fontsize <- set_fontsize(mat)

    h <- Heatmap(mat,
        column_title = "Promoter mean methylation in embryos",
        name = "Meth",
        col = colorRamp2(seq_len(100),
                         colorRampPalette(c("moccasin", "dodgerblue4"))(100)),
        na_col = "gray80",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_heatmap_legend = TRUE,
        show_row_dend = FALSE,
        column_split = df_col$Stage,
        column_gap = unit(c(1), "mm"),
        top_annotation = c(column_ha_stage,
                           column_ha_genotype,
          gap = unit(0.1, "mm")), 
        row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = 0),
        column_names_side = "bottom",
        row_names_side = "right")

    ifelse(values_only, return(database[, df_col$Sample_Name]), return(h))
}
