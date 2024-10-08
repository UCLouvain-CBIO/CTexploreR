#' Methylation of CpGs located in promoters in normal tissues
#'
#' @description Plots a heatmap of the methylation of CpGs located in a
#'  promoter, in normal tissues. X-axis corresponds to the
#' CpGs position (related to TSS).
#'
#' @param gene Name of selected gene
#'
#' @param nt_up Number of nucleotides upstream the TSS to analyse
#' (by default 1000, maximum value 5000)
#'
#' @param nt_down Number of nucleotides downstream the TSS to analyse
#' (by default 200, maximum value 5000)
#'
#' @param values_only Boolean (FALSE by default). If set to TRUE, the function 
#' will return the methylation values of all cytosines in the promoter instead 
#' of the heatmap.
#'
#' @return Heatmap of the methylation of CpGs located in a 
#' promoter, in normal tissues. If `values_only` = TRUE, methylation values are
#' returned instead.
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
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#'
#' @examples
#' normal_tissues_methylation(gene = "TDRD1", 1000, 0)
normal_tissues_methylation <- function(gene, nt_up = 1000, nt_down = 200,
                                       values_only = FALSE) {
    suppressMessages({
        database <- CTdata::methylation_in_tissues()
        all_genes <- CTdata::all_genes()
    })

    if (!gene %in% all_genes$external_gene_name) {
        stop(gene, " is not in the CT database")
    }

    chr <- all_genes |>
        filter(.data$external_gene_name == gene) |>
        pull(.data$chr)

    strand <- all_genes |>
        filter(.data$external_gene_name == gene) |>
        pull(.data$strand)

    TSS <- all_genes |>
        filter(.data$external_gene_name == gene) |>
        pull(.data$transcription_start_site)
    
    if (nt_up > 5000) {
      nt_up <- 5000
      warning("Replacing `nt_up` value by 5000, the maximum number of 
      nucleotides upstream the TSS analysable")
    }
    
    if (nt_down > 5000) {
      nt_down <- 5000
      warning("Replacing `nt_down` value by 5000, the maximum number of 
      nucleotides downstream the TSS analysable")
    }

    if (strand == 1) {
        promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                               strand = "+",
                               ranges = IRanges(
                                 start = TSS - nt_up,
                                 end = TSS + nt_down))
        promoter_gr$TSS <- TSS
    }

    if (strand == -1) {
        promoter_gr <- GRanges(seqnames = paste0("chr", chr),
                               strand = "-",
                               ranges = IRanges(start = TSS - nt_down,
                                                end = TSS + nt_up))
        promoter_gr$TSS <- TSS
    }

    promoter_methylation <- subsetByOverlaps(database, promoter_gr)

    methylation_individual_CpG <- suppressWarnings(
        as_tibble(assay(promoter_methylation)) |>
            mutate(TSS = promoter_gr$TSS) |>
            mutate(CG_pos = promoter_methylation@rowRanges@ranges@start) |>
            mutate(relative_pos = case_when(
                strand == 1 ~ .data$CG_pos - TSS,
                strand == -1 ~ TSS - .data$CG_pos)) |>
            dplyr::select(.data$CG_pos, .data$relative_pos,
                colnames(promoter_methylation)) |>
            pivot_longer(!c(.data$CG_pos, .data$relative_pos),
                         names_to = "cell",
                         values_to = "met") |>
            dplyr::select(-.data$CG_pos) |>
            pivot_wider(names_from = .data$relative_pos, 
                        values_from = .data$met))

    mat <- as.matrix(methylation_individual_CpG[, -1])
    rownames(mat) <- methylation_individual_CpG$cell

    h <- Heatmap(mat,
        name = "Meth",
        col = colorRamp2(seq_len(100),
                         colorRampPalette(c("moccasin", "dodgerblue4"))(100)),
        na_col = "gray80",
        column_title = paste0(gene, " (TSS chr", chr, ":", TSS, ")"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_heatmap_legend = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        row_names_side = "left")

    ifelse(values_only, return(mat), return(h))
}
