#' Methylation-Expression correlation of Cancer-Testis genes in TCGA samples
#'
#' @description Plots the correlation between methylation and
#' expression values of a Cancer-Testis (CT) gene in TCGA samples.
#'
#' @param gene CT gene name(s)
#'
#' @param tumor TCGA tumor code. c("SKCM", "LUAD", "LUSC", "COAD", "ESCA",
#' "BRCA", "HNSC", "all")
#'
#' @param corr_coeff Boolean (FALSE by default). If set to TRUE, the function
#' will invisibly return the correlation coefficient (Pearson), between
#' methylation and expression values for the gene in selected samples.
#'
#' @param return Boolean (FALSE by default). If set to TRUE, the function will
#' return the methylation and expression values in all samples instead of the
#' heatmap.
#'
#' @details The coefficient of correlation is set to `NA` if no probes are
#' found in promoter regions or if less than 1% of tumors are positive
#' (TPM >= 1) for the gene.
#'
#' @return A correlation plot between gene expression and methylation
#' values of probe(s) located in its promoter region (defined as 1000
#' nucleotides upstream TSS and 200 nucleotides downstream TSS).
#' If return = TRUE, methylation and expression values for the gene in
#' selected tumors are returned in a tibble instead. If `corr_coeff` is
#' set to TRUE, the correlation coefficient is being returned instead.
#'
#' @export
#'
#' @importFrom BiocGenerics intersect
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#' @importFrom tibble tibble enframe
#' @importFrom dplyr left_join select mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth ggtitle xlim
#' @importFrom stats cor.test
#'
#' @examples
#' TCGA_methylation_expression_correlation(tumor = "LUAD", gene = "TDRD1")
#' TCGA_methylation_expression_correlation(tumor = "all", gene = "TDRD1")
#' TCGA_methylation_expression_correlation(tumor = "all", gene = "TDRD1",
#' corr_coeff = TRUE)
# TCGA_methylation_expression_correlation <- function(tumor, gene = NULL,
#                                                     return = FALSE,
#                                                     corr_coeff = FALSE) {
# 
#   if (missing(tumor)) {
#     stop("TCGA tumor code ('SKCM', 'LUAD', 'LUSC', 'COAD', 'ESCA', 'BRCA', 'HNSC' or 'all') must be specified!")
#   }
# 
#   if (!all(tumor %in% c("SKCM", "LUAD", "LUSC", "COAD", "ESCA", "BRCA",
#                         "HNSC", "all"))) {
#     stop("TCGA tumor code must be one of ('SKCM', 'LUAD', 'LUSC', 'COAD', 'ESCA', 'BRCA', 'HNSC', or 'all')!")
#   }
# 
#   if (is.null(gene)) {
#     stop("Gene name must be specified!")
#   }
# 
#   if (length(tumor) == 1 & tumor[1] == "all") {
#     tumor <- c('SKCM', 'LUAD', 'LUSC', 'COAD', 'ESCA', 'BRCA', 'HNSC')
#   }
# 
#     TPM <- TCGA_TPM
#     met <- TCGA_CT_methylation
# 
#   if (!gene %in% CT_genes$external_gene_name) { ## Warning for length(gene) > 1
#     message(paste0("Check gene name! ", gene, " is not in the database"))
#     return(invisible(NA))
#   }
# 
#   ## Load expression and methylation data from selected tumors
#   ##(keeping only primary and metastatic tumors)
#   projects <- paste0("TCGA-", tumor)
#   TPM <- TPM[, colData(TPM)$project_id %in% projects]
#   TPM <- TPM[, colData(TPM)$shortLetterCode == "TP" |
#                colData(TPM)$shortLetterCode == "TM" |
#                colData(TPM)$shortLetterCode == "NT"]
# 
#   met <- met[, colData(met)$project_id %in% projects]
# 
#   ## select tumors for which expression and methylation data are available
#   samples <- intersect(colData(TPM)$sample, colData(met)$sample)
# 
#   ## Create a Grange corresponding to promoter regions
#   nt_up <- 1000
#   nt_down <- 200
#   CT_promoter_gr <- makeGRangesFromDataFrame(
#     CT_genes %>%
#       dplyr::select(ensembl_gene_id, external_gene_name,
#                     external_transcript_name,chromosome_name, strand,
#                     transcription_start_site) %>%
#       mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
#       mutate(strand = ifelse(strand == 1, '+', '-')) %>%
#       mutate(start = case_when(strand == '+' ~ transcription_start_site - nt_up,
#                                strand == '-' ~ transcription_start_site - nt_down)) %>%
#       mutate(stop = case_when(strand == '+' ~ transcription_start_site + nt_down,
#                               strand == '-' ~ transcription_start_site + nt_up)),
#     keep.extra.columns = TRUE,
#     seqnames.field = "chromosome_name",
#     start.field = "start",
#     end.field = "stop")
# 
#   met_roi <- subsetByOverlaps(met[, colData(met)$sample %in% samples],
#                               CT_promoter_gr[CT_promoter_gr$external_gene_name == gene])
# 
#   ## Evaluate mean methylation value of promoter probe(s) in each sample
#   met_mean <- colMeans(assay(met_roi), na.rm = TRUE)
#   names(met_mean) <- substr(names(met_mean), 1, 16)
# 
#   ensembl <- CT_genes[CT_genes$external_gene_name == gene, "ensembl_gene_id",
#                       drop = TRUE]
#   TPM <- assay(TPM[rownames(TPM) %in% ensembl,
#                    colData(TPM)$sample %in% samples])
# 
#   colnames(TPM) <- substr(colnames(TPM), 1, 16)
# 
#   if (nrow(TPM) == 0) {
#     print(paste0(gene, " is not in TCGA expression database"))
#     return(invisible(NA))
#   }
# 
#   suppressMessages(
#     methylation_expression <- left_join(enframe(met_mean, name = "sample",
#                                                 value = "met"),
#                                         enframe(TPM[1,], name = "sample",
#                                                 value = "TPM"))
#   )
# 
#   suppressMessages(
#     methylation_expression <- methylation_expression %>%
#     left_join(as_tibble(colData(met)) %>%
#                 mutate(Tumor = sub(pattern = "TCGA-", x = project_id,
#                                    replacement = '')) %>%
#                 dplyr::select(sample, Tumor))
#   )
# 
#   ## stop if no probes or no methylation values for probes within the region
#   if (all(is.na(methylation_expression$met))) {
#     message(paste0("No probes for ", gene))
#     return(invisible(NA))
#   }
# 
#   ## Gene has to be expressed (TPM >= 1) in at least 1%
#   ## of the samples to evaluate correlation
#   if (length(methylation_expression$TPM[methylation_expression$TPM] >= 1) <
#       ceiling(nrow(methylation_expression) * 0.01)) {
#     message(paste0("Too few positive samples to estimate a correlation for ",
#                    gene))
#     cor <- NA
#   } else {
#     cor <- cor.test(methylation_expression$met,
#                     log1p(methylation_expression$TPM))$estimate
#   }
# 
#   TPM <- met <- NULL
# 
#   methylation_expression$Tissue <- ifelse(substr(methylation_expression$sample,
#                                                  14, 15) == "11",
#                                           "Peritumoral", "Tumor")
# 
#   if (length(unique(methylation_expression$Tumor)) > 1) {  # color by tumor type
#     p <- ggplot(methylation_expression[order(methylation_expression$Tissue,
#                                              decreasing = TRUE), ],
#                 aes(x = met, y = log1p(TPM))) +
#       geom_point(alpha = 0.6, aes(color = Tumor, shape = Tissue)) +
#       ggtitle(paste0(gene, "(Pearson's corr = ", round(cor, 2), ")")) +
#       xlim(0, 1)
#   }
# 
#   if (length(unique(methylation_expression$Tumor)) == 1) {  # color by tissue if only one tumor type
#     p <- ggplot(methylation_expression[order(methylation_expression$Tissue,
#                                              decreasing = TRUE), ],
#                 aes(x = met, y = log1p(TPM))) +
#       geom_point(alpha = 0.6, aes(color = Tissue)) +
#       ggtitle(paste0(gene, " in ", tumor, " (corr = ", round(cor, 2), ")")) +
#       xlim(0, 1)
#   }
# 
# 
#   if (corr_coeff) {
#     unname(round(cor, 2))
#   } else {
#     if (return == FALSE) {
#       suppressMessages(print(p))
#     } else {
#       methylation_expression
#     }
# 
#   }
# 
# }


TCGA_methylation_expression_correlation <- function(tumor, gene = NULL,
                                                    return = FALSE,
                                                    corr_coeff = FALSE) {
  
  TPM <- TCGA_TPM
  TPM$type <- sub(pattern = "TCGA-", x = colData(TPM)$project_id, replacement = '')
  met <- TCGA_CT_methylation
  met$type <- sub(pattern = "TCGA-", x = colData(met)$project_id, replacement = '')
  
  valid_tumor <- c(unique(colData(TPM)$type), "all")
  type <- check_names(variable = tumor, valid_vector = valid_tumor)
  stopifnot("No valid tumor type entered" = length(type) > 0)
  if (type != "all") {
    TPM <- TPM[, TPM$type %in% type]
    met <- met[, met$type %in% type]
  } 
  
  stopifnot("No valid gene name entered" = !is.null(gene))
  valid_gene_names <- CT_genes$external_gene_name
  gene <- check_names(gene, valid_gene_names)
  TPM <- TPM[rowData(TPM)$external_gene_name %in% gene, ]
  
  
  
  
  
  # if (!gene %in% CT_genes$external_gene_name) { ## Warning for length(gene) > 1
  #   message(paste0("Check gene name! ", gene, " is not in the database"))
  #   return(invisible(NA))
  # }
  
  ## select tumors for which expression and methylation data are available
  samples <- intersect(colData(TPM)$sample, colData(met)$sample)
  
  ## Create a Grange corresponding to promoter regions
  nt_up <- 1000
  nt_down <- 200
  CT_promoter_gr <- makeGRangesFromDataFrame(
    CT_genes %>%
      dplyr::select(ensembl_gene_id, external_gene_name,
                    external_transcript_name,chromosome_name, strand,
                    transcription_start_site) %>%
      mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
      mutate(strand = ifelse(strand == 1, '+', '-')) %>%
      mutate(start = case_when(strand == '+' ~ transcription_start_site - nt_up,
                               strand == '-' ~ transcription_start_site - nt_down)) %>%
      mutate(stop = case_when(strand == '+' ~ transcription_start_site + nt_down,
                              strand == '-' ~ transcription_start_site + nt_up)),
    keep.extra.columns = TRUE,
    seqnames.field = "chromosome_name",
    start.field = "start",
    end.field = "stop")
  
  met_roi <- subsetByOverlaps(met[, colData(met)$sample %in% samples],
                              CT_promoter_gr[CT_promoter_gr$external_gene_name == gene])
  
  ## Evaluate mean methylation value of promoter probe(s) in each sample
  met_mean <- colMeans(assay(met_roi), na.rm = TRUE)
  names(met_mean) <- substr(names(met_mean), 1, 16)
  
  ensembl <- CT_genes[CT_genes$external_gene_name == gene, "ensembl_gene_id",
                      drop = TRUE]
  TPM <- assay(TPM[rownames(TPM) %in% ensembl,
                   colData(TPM)$sample %in% samples])
  
  colnames(TPM) <- substr(colnames(TPM), 1, 16)
  
  if (nrow(TPM) == 0) {
    print(paste0(gene, " is not in TCGA expression database"))
    return(invisible(NA))
  }
  
  suppressMessages(
    methylation_expression <- left_join(enframe(met_mean, name = "sample",
                                                value = "met"),
                                        enframe(TPM[1,], name = "sample",
                                                value = "TPM"))
  )
  
  suppressMessages(
    methylation_expression <- methylation_expression %>%
      left_join(as_tibble(colData(met)) %>%
                  mutate(Tumor = sub(pattern = "TCGA-", x = project_id,
                                     replacement = '')) %>%
                  dplyr::select(sample, Tumor))
  )
  
  ## stop if no probes or no methylation values for probes within the region
  if (all(is.na(methylation_expression$met))) {
    message(paste0("No probes for ", gene))
    return(invisible(NA))
  }
  
  ## Gene has to be expressed (TPM >= 1) in at least 1%
  ## of the samples to evaluate correlation
  if (length(methylation_expression$TPM[methylation_expression$TPM] >= 1) <
      ceiling(nrow(methylation_expression) * 0.01)) {
    message(paste0("Too few positive samples to estimate a correlation for ",
                   gene))
    cor <- NA
  } else {
    cor <- cor.test(methylation_expression$met,
                    log1p(methylation_expression$TPM))$estimate
  }
  
  TPM <- met <- NULL
  
  methylation_expression$Tissue <- ifelse(substr(methylation_expression$sample,
                                                 14, 15) == "11",
                                          "Peritumoral", "Tumor")
  
  if (length(unique(methylation_expression$Tumor)) > 1) {  # color by tumor type
    p <- ggplot(methylation_expression[order(methylation_expression$Tissue,
                                             decreasing = TRUE), ],
                aes(x = met, y = log1p(TPM))) +
      geom_point(alpha = 0.6, aes(color = Tumor, shape = Tissue)) +
      ggtitle(paste0(gene, "(Pearson's corr = ", round(cor, 2), ")")) +
      xlim(0, 1)
  }
  
  if (length(unique(methylation_expression$Tumor)) == 1) {  # color by tissue if only one tumor type
    p <- ggplot(methylation_expression[order(methylation_expression$Tissue,
                                             decreasing = TRUE), ],
                aes(x = met, y = log1p(TPM))) +
      geom_point(alpha = 0.6, aes(color = Tissue)) +
      ggtitle(paste0(gene, " in ", tumor, " (corr = ", round(cor, 2), ")")) +
      xlim(0, 1)
  }
  
  
  if (corr_coeff) {
    unname(round(cor, 2))
  } else {
    if (return == FALSE) {
      suppressMessages(print(p))
    } else {
      methylation_expression
    }
    
  }
  
}



