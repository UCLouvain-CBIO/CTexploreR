#' Promoter methylation of Cancer-Testis genes in normal tissues
#'
#' @description Plots a heatmap of Promoter methylation of a Cancer-Testis (CT)
#' gene in normal tissues
#'
#' @param database CT_methylation_in_tissues
#'
#' @param gene name of selected CT gene
#'
#' @param nt_up number of nucleotides upstream the TSS to analyse (1000 by default)
#'
#' @param nt_down number of nucleotides downstream the TSS to analyse (200 by default)
#'
#' @return heatmap of Promoter methylation of a Cancer-Testis (CT) gene in normal tissues.
#' Methylation values are returned invisibly.
#'
#' @export
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom tibble as_tibble
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom dplyr filter pull mutate select case_when
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom magrittr %>%
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' normal_tissues_methylation(CT_methylation_in_tissues, gene = "TDRD1", 1000, 0)
normal_tissues_methylation <- function(database, gene, nt_up = NULL, nt_down = NULL){

  if (missing(database)) {
    stop("Database must be specified!")
  }

  if (!missing(database)) {
    data <- database
  }

  if (!gene %in% CT_genes$external_gene_name) {
    stop(paste0(gene, " is not in the CT database"))
  }

  chr <- CT_genes %>%
    filter(external_gene_name == gene) %>%
    pull(chromosome_name)

  strand <- CT_genes %>%
    filter(external_gene_name == gene) %>%
    pull(strand)

  TSS <- CT_genes %>%
    filter(external_gene_name == gene) %>%
    pull(transcription_start_site)

  if (is.null(nt_up)) {nt_up <- 1000}

  if (is.null(nt_down)) {nt_down <- 200}

  if (strand == 1) {
    promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                           strand = '+',
                           ranges = IRanges(start = TSS - nt_up, end = TSS + nt_down))
    promoter_gr$TSS <- TSS
  }

  if (strand == -1) {
    promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                           strand = '-',
                           ranges = IRanges(start = TSS - nt_down, end = TSS + nt_up))
    promoter_gr$TSS <- TSS
  }

  promoter_methylation <- subsetByOverlaps(data, promoter_gr)

  methylation_individual_CpG <- suppressWarnings(
    as_tibble(assay(promoter_methylation)) %>%
      mutate(TSS = promoter_gr$TSS) %>%
      mutate(CG_pos = promoter_methylation@rowRanges@ranges@start) %>%
      mutate(relative_pos = case_when(strand == 1 ~ CG_pos - TSS,
                                      strand == -1 ~ TSS - CG_pos)) %>%
      dplyr::select(CG_pos, relative_pos, colnames(promoter_methylation)) %>%
      pivot_longer(!c(CG_pos, relative_pos), names_to = "cell", values_to = "met") %>%
      dplyr::select(-CG_pos) %>%
      pivot_wider(names_from = relative_pos, values_from = met))

  mat <- as.matrix(methylation_individual_CpG[,-1])
  rownames(mat) <- methylation_individual_CpG$cell

  h <- Heatmap(mat,
               name=paste0('methylation'),
               col = colorRamp2(c(1:100),
                                colorRampPalette(c("moccasin","dodgerblue4"))(100)),
               na_col = "gray80",
               column_title = paste0(gene, " (TSS chr", chr, ":", TSS, ")"),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               show_heatmap_legend = TRUE,
               row_names_gp = gpar(fontsize = 8),
               column_names_side = "bottom",
               row_names_side = "left")

  print(h)
  invisible(mat)
}

