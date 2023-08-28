#' Expression values (TPM) of genes in normal tissues with or without multimapping
#'
#' @description
#'
#' Plots a heatmap of gene expression values in a set of normal
#' tissues. Expression values (in TPM) have been evaluated by either
#' counting or discarding multi-mapped reads. Indeed, many CT genes
#' belong to gene families from which members have identical or nearly
#' identical sequences. Some CT can only be detected in RNAseq data in
#' which multimapping reads are not discarded.
#'
#' @param genes `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT genes.
#'
#' @param multimapping `logical(1)` that specifies if returned
#'     expression values must take into account or not multi-mapped
#'     reads
#'
#' @param units `character(1)` with expression values unit.  Can be
#'     `"TPM"` (default) or `"log_TPM"` (log(TPM + 1)).
#'
#' @param return `logical(1)`. If `TRUE`, the function will return the
#'     expression values in all samples instead of the
#'     heatmap. Default is `FALSE`.
#'
#' @details
#'
#' RNAseq data from a set of normal tissues were downloaded from
#' Encode.  (see inst/scripts/make_CT_normal_tissues_multimapping.R
#' for fastq references) Fastq files were processed using a standard
#' RNAseq pipeline including
#' [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#' for the quality control of the raw data, and
#' [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to
#' remove low quality reads and trim the adapter from the sequences.
#' [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used
#' to align reads to grch38 genome.
#' [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html)
#' was used to assign reads to genes using
#' Homo_sapiens.GRCh38.105.gtf.
#'
#' Two different pipelines were run in order to remove or not
#' multi-mapping reads.  When multimapping was allowed, hisat2 was run
#' with -k 20 parameter (reports up to 20 alignments per read), and
#' featurecounts was run with -M parameter (multi-mapping reads are
#' counted).
#'
#' @return A heatmap of selected gene expression values in a set of
#'     normal tissues calculated by counting or discarding
#'     multi-mapped reads.  If `return = TRUE`, gene expression values
#'     are returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' normal_tissue_expression_multimapping(
#'    genes = c("GAGE13", "CT45A6", "NXF2","SSX2", "CTAG1A", "MAGEA3", "MAGEA6"),
#'    multimapping = FALSE)
#' normal_tissue_expression_multimapping(
#'    genes = c("GAGE13", "CT45A6", "NXF2", "SSX2", "CTAG1A", "MAGEA3", "MAGEA6"),
#'    multimapping = TRUE)
normal_tissue_expression_multimapping <- function(genes = NULL,
                                                  multimapping = NULL,
                                                  units = "TPM",
                                                  return = FALSE) {
    if (is.null(multimapping))
        stop("multimapping parameter should be set to TRUE/FALSE")

    suppressMessages({
      CT_genes <- CTdata::CT_genes()
      database <- CTdata::normal_tissues_multimapping_data()
    })

    if (is.null(genes)) genes <- CT_genes$external_gene_name
    valid_gene_names <- unique(rowData(database)$external_gene_name)
    genes <- check_names(genes, valid_gene_names)
    database <- database[rowData(database)$external_gene_name %in% genes, ]

    if (multimapping) {
        mat <- assay(database, "TPM_with_multimapping")
        title <- "Expression (multi-mapped reads were counted)"
    } else {
        mat <- assay(database, "TPM_no_multimapping")
        title <- "Expression (multimapped reads were discared)"
    }
    rownames(mat) <- rowData(database)$external_gene_name

    name <- "TPM"
    if (units == "log_TPM") {
        mat <- log1p(mat)
        name <- "log_TPM"
    }

    if (dim(mat)[1] > 100) fontsize <- 4
    if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) fontsize <- 5
    if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) fontsize <- 6
    if (dim(mat)[1] <= 20) fontsize <- 8

    h <- Heatmap(mat,
                 name = name,
                 col = colorRamp2(seq(0, max(mat), length = 11),
                                  legend_colors),
                 column_title = title,
                 cluster_rows = TRUE,
                 show_row_dend = FALSE,
                 cluster_columns = TRUE,
                 show_column_dend = FALSE,
                 row_names_gp = gpar(fontsize = fontsize),
                 column_names_gp = gpar(fontsize = 6),
                 clustering_method_rows = "ward.D")


    if (return)
        return(mat)

    return(h)
}
