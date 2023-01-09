#' Expression values (TPM) of genes in normal tissues with or without multimapping
#'
#' @description Plots a heatmap of gene expression values in a set of normal
#' tissues. Expression values (in TPM) have been evaluated by either counting or
#' discarding multi-mapped reads. Indeed, many CT genes belong to gene families
#' from which members have identical or nearly identical sequences. Some CT can
#' only be detected in RNAseq data in which multimapping reads are not discared.
#'
#' @param database normal_tissues_multimapping_data
#'
#' @param genes Genes selected (all CT genes by default)
#'
#' @param multimapping Set to TRUE or FALSE to specify if returned expression
#' values must take into account or not multi-mapped reads
#'
#' @param units Expression values units.
#' Can be "TPM" (default) or "log_TPM" (log(TPM + 1))
#'
#' @details
#' RNAseq data from a set of normal tissues were downloaded from Encode.
#' (see inst/scripts/make_CT_normal_tissues_multimapping.R for fastq references)
#' Fastq files were processed using a standard RNAseq pipeline including
#' [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for the
#' quality control of the raw data, and
#' [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
#' to remove low quality reads and trim the adapter from the sequences.
#' [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) was used to align
#' reads to grch38 genome.
#' [featurecounts](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) was used
#' to assign reads to genes using Homo_sapiens.GRCh38.105.gtf.
#'
#' Two different pipelines were run in order to remove or not multi-mapping reads.
#' When multimapping was allowed, hisat2 was run with -k 20 parameter (reports
#' up to 20 alignments per read), and featurecounts was run with -M parameter
#' (multi-mapping reads are counted).
#'
#' @return A heatmap of selected gene expression values (TPM) in a set of
#' normal tissues calculated by counting or discarding multi-mapped reads.
#' Gene TPM values are invisibly returned.
#'
#' @export
#'
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
#' @examples
#' normal_tissue_expression_multimapping(database = normal_tissues_multimapping_data,
#' genes = c("GAGE13", "CT45A6", "NXF2", "SSX2", "CTAG1A", "MAGEA3", "MAGEA6"),
#' multimapping = FALSE)
#' normal_tissue_expression_multimapping(database = normal_tissues_multimapping_data,
#' genes = c("GAGE13", "CT45A6", "NXF2", "SSX2", "CTAG1A", "MAGEA3", "MAGEA6"),
#' multimapping = TRUE)
normal_tissue_expression_multimapping <- function(database, genes = NULL,
                                                  multimapping = NULL,
                                                  units = "TPM") {

  if (missing(database)) {
    stop("Database must be specified!")
  }

  if (!missing(database)) {
    data <- database
  }

  if (is.null(multimapping)) {
    stop("multimapping parameter should be set to TRUE/FALSE")
  }

  if (multimapping == TRUE) {
    mat <- assay(data, "TPM_with_multimapping")
    title <- "Gene expression in normal tissues (multimapped reads were counted)"
  }

  if (multimapping == FALSE) {
    mat <- assay(data, "TPM_no_multimapping")
    title <- "Gene expression in normal tissues (multimapped reads were discared)"
  }
  rownames(mat) <- rowData(data)$external_gene_name

  if (is.null(genes)) {
    mat <- mat[CT_genes$external_gene_name, ]
  }

  if (!is.null(genes)) {
    if (!all(genes %in% rownames(mat))) {
      message("Check gene name(s)!\n")
      message(paste0(genes[!genes %in% rownames(mat)],
                     " is not in the database.\n"))
    }
    mat <- mat[genes[genes %in% rownames(mat)], , drop = FALSE]
  }

  name <- "TPM"
  if (units == "log_TPM") {
    mat <- log1p(mat)
    name <- "log_TPM"
  }

  if (dim(mat)[1] > 100) { fontsize <- 4 }
  if (dim(mat)[1] > 50 & dim(mat)[1] <= 100) { fontsize <- 5 }
  if (dim(mat)[1] > 20 & dim(mat)[1] <= 50) { fontsize <- 6 }
  if (dim(mat)[1] <= 20) { fontsize <- 8 }

  h <- suppressMessages(Heatmap(mat,
                                name = name,
                                col = colorRamp2(seq(0, max(mat), length = 11),
                                                 c("#5E4FA2", "#3288BD",
                                                   "#66C2A5", "#ABDDA4",
                                                   "#E6F598", "#FFFFBF",
                                                   "#FEE08B", "#FDAE61",
                                                   "#F46D43", "#D53E4F",
                                                   "#9E0142")),
                                column_title = title,
                                cluster_rows = TRUE,
                                show_row_dend = FALSE,
                                cluster_columns = TRUE,
                                show_column_dend = FALSE,
                                row_names_gp = gpar(fontsize = fontsize),
                                column_names_gp = gpar(fontsize = 6),
                                clustering_method_rows = "ward.D"))

  print(h)
  invisible(mat)
}
