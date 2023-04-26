#' Methylation-Expression correlation of Cancer-Testis genes in TCGA samples
#'
#' @description
#'
#' Plots the correlation between methylation and expression values of
#' a Cancer-Testis (CT) gene in TCGA samples.
#'
#' @param gene `character` nameing the selected genes. The default
#'     value, `NULL`, takes all CT genes.
#'
#' @param tumor `character` defining the TCGA tumor type. Can be one
#'     of "SKCM", "LUAD", "LUSC", "COAD", "ESCA", "BRCA", "HNSC", or
#'     "all" (default).
#'
#' @param nt_up `numeric(1)` indicating the number of nucleotides
#'     upstream the TSS to define the promoter region (1000 by
#'     default)
#'
#' @param nt_down `numeric(1)` indicating the number of nucleotides
#'     downstream the TSS to define the promoter region (200 by
#'     default)
#'
#' @param corr_coeff `logical(1)`. If `TRUE`, the function will
#'     invisibly return the correlation coefficient (Pearson), between
#'     methylation and expression values for the gene in selected
#'     samples. Default is `FALSE`.
#'
#' @param return `logical(1)`. If `TRUE`, the function will return the
#'     methylation and expression values in all samples instead of the
#'     heatmap. Default is `FALSE`.
#'
#' @details The coefficient of correlation is set to `NA` if no probes
#'     are found in promoter regions or if less than 1% of tumors are
#'     positive (TPM >= 1) for the gene.
#'
#' @return A correlation plot between gene expression and methylation
#'     values of probe(s) located in its promoter region (defined as
#'     1000 nucleotides upstream TSS and 200 nucleotides downstream
#'     TSS).  If return = TRUE, methylation and expression values for
#'     the gene in selected tumors are returned in a tibble
#'     instead. If `corr_coeff` is set to TRUE, the correlation
#'     coefficient is being returned instead.
#'
#' @export
#'
#' @importFrom BiocGenerics intersect
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#' @importFrom tibble tibble
#' @importFrom dplyr filter left_join select mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth ggtitle xlim
#' @importFrom stats cor.test quantile
#' @importFrom rlang .data
#'
#' @examples
#' TCGA_methylation_expression_correlation("LUAD", gene = "TDRD1",
#'                                         return = FALSE)
#' TCGA_methylation_expression_correlation(c("LUAD", "LUSC"),
#'                                         gene = "MAGEA1", return = FALSE)
#' TCGA_methylation_expression_correlation("LUAD",
#'                                         gene = "TDRD1", return = TRUE)
TCGA_methylation_expression_correlation <- function(tumor,
                                                    gene = NULL,
                                                    nt_up = 1000,
                                                    nt_down = 200,
                                                    return = FALSE,
                                                    corr_coeff = FALSE) {
    suppressMessages({
        CT_genes <- CTdata::CT_genes()
        TPM <- CTdata::TCGA_TPM()
        met <- CTdata::TCGA_CT_methylation()
    })

    TPM$type <- sub(pattern = "TCGA-", x = colData(TPM)$project_id,
                    replacement = '')
    met$type <- sub(pattern = "TCGA-", x = colData(met)$project_id,
                    replacement = '')

    valid_tumor <- c(unique(colData(TPM)$type), "all")
    type <- check_names(variable = tumor, valid_vector = valid_tumor)
    stopifnot("No valid tumor type entered" = length(type) > 0)
    if (!"all" %in% type) {
        TPM <- TPM[, TPM$type %in% type]
        met <- met[, met$type %in% type]
    }

    stopifnot("No valid gene name entered" = !is.null(gene))
    valid_gene_names <- CT_genes$external_gene_name
    valid_gene_names <- valid_gene_names[valid_gene_names %in% rowData(TPM)$external_gene_name]
    gene <- check_names(gene, valid_gene_names)
    stopifnot("Gene not in TCGA database" = gene %in% rowData(TPM)$external_gene_name)
    stopifnot(length(gene) == 1)
    TPM <- TPM[rowData(TPM)$external_gene_name %in% gene, ]

    ## select tumors for which expression and methylation data are available
    samples <- intersect(colData(TPM)$sample, colData(met)$sample)
    TPM <- TPM[, TPM$sample %in% samples]
    met <- met[, met$sample %in% samples]

    ## Create a Grange corresponding to promoter region (defined as nt_up
    ## nucleotides upstream and nt_down nucleotides downstream TSS)
    ## Calculates mean methylation value of promoter probe(s) in each sample
    gene_promoter_gr <- makeGRangesFromDataFrame(
        CT_genes |>
        dplyr::filter(.data$external_gene_name == gene) |>
        dplyr::select("ensembl_gene_id", "external_gene_name",
                      "external_transcript_name", "chr", "strand",
                      "transcription_start_site") |>
        mutate(chr = paste0("chr", .data$chr)) |>
        mutate(strand = ifelse(.data$strand == 1, '+', '-')) |>
        mutate(start = case_when(strand == '+' ~ .data$transcription_start_site - nt_up,
                                 strand == '-' ~ .data$transcription_start_site - nt_down)) |>
        mutate(stop = case_when(strand == '+' ~ .data$transcription_start_site + nt_down,
                                strand == '-' ~ .data$transcription_start_site + nt_up)),
        keep.extra.columns = TRUE,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "stop")

    met_roi <- subsetByOverlaps(met, gene_promoter_gr)
    met_mean <- colMeans(assay(met_roi), na.rm = TRUE)

    ## Keep prefix of TCGA sample names (TCGA-XX-XXXX-XXX) to be able to join
    ## expression and methylation data
    names(met_mean) <- substr(names(met_mean), 1, 16)
    colnames(TPM) <- substr(colnames(TPM), 1, 16)

    ## Rm duplicated samples
    TPM <- TPM[, !duplicated(colnames(TPM))]
    met_mean <- met_mean[!duplicated(names(met_mean))]
    
    methylation_expression <-
      suppressMessages(left_join(
        tibble(sample = names(met_mean), met = met_mean),
        tibble(sample = colnames(TPM), TPM = as.vector(assay(TPM)))))

    ## stop if no probes or no methylation values for probes within the region
    if (all(is.na(methylation_expression$met))) {
        message(paste0("No probes or no methylation values for probes for ", gene))
        cor <- NA
    }

    ## Gene has to be expressed (TPM >= 1) in at least 1%
    ## of the samples to evaluate correlation
    if (quantile(methylation_expression$TPM, 0.99) < 1) {
        message(paste0("Too few positive samples to estimate a correlation for ",
                       gene))
        cor <- NA
    } else {
        cor <- cor.test(methylation_expression$met,
                        log1p(methylation_expression$TPM))$estimate
    }

    methylation_expression <-
        merge(methylation_expression,
              colData(TPM)[, c("sample", "shortLetterCode", "type")])

    ## shortLetterCode TP and TM <=> primary / metastatic tumors
    ## shortLetterCode NT <=> Normal peritumoral tissue
    methylation_expression$Tissue <- "Tumor"
    methylation_expression$Tissue[methylation_expression$shortLetterCode == "NT"] <- "Peritumoral"

    ## Color by tissue (Peritumoral / Tumor) if only one tumor type,
    ## otherwise color by tumor types.
    if (all(type == "all") | length(type) > 1) {
        p <- ggplot(
            as_tibble(methylation_expression[order(methylation_expression$type,
                                                   decreasing = TRUE), ]),
            aes(x = .data$met, y = log1p(.data$TPM))) +
            geom_point(alpha = 0.8, aes(color = .data$type, shape = .data$Tissue)) +
            ggtitle(paste0(gene, "(Pearson's corr = ", round(cor, 2), ")")) +
            scale_colour_manual(values = TCGA_colors) +
            xlim(0, 1)
    } else {
        p <- ggplot(
            as_tibble(methylation_expression[order(methylation_expression$Tissue,
                                                   decreasing = TRUE), ]),
            aes(x = .data$met, y = log1p(.data$TPM))) +
            geom_point(alpha = 0.8, aes(color = .data$Tissue)) +
            ggtitle(paste0(gene, " in ", tumor, " (corr = ", round(cor, 2), ")")) +
            scale_colour_manual(values = TCGA_colors) +
            xlim(0, 1)
    }

    if (corr_coeff)
        return(unname(cor))
    if (return)
        return(methylation_expression)

    return(p)
}
