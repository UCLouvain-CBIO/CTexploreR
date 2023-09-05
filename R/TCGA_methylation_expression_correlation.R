#' Methylation-Expression correlation of Cancer-Testis genes in TCGA
#' samples
#'
#' @description
#'
#' Plots the correlation between methylation and expression values of
#' a Cancer-Testis (CT) gene in TCGA samples.
#'
#' @param gene `character` selected gene. 
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
#' @param min_probe_number `numeric(1)` indicating the minimum 
#'     number of probes (with methylation values) within the selected region 
#'     to calculate its mean methylation level. Default is 3.
#'     
#' @param include_normal_tissues `logical(1)`. If `TRUE`, the function will 
#'     include normal peritumoral tissues in addition to tumoral samples. 
#'     Default is `FALSE`.
#'
#' @param return `logical(1)`. If `TRUE`, the function will return the
#'     methylation and expression values in TCGA samples instead of the
#'     heatmap. Default is `FALSE`.
#'
#' @details The coefficient of correlation is set to `NA` if no probes
#'     are found in promoter regions or if less than 1% of tumors are
#'     positive (TPM >= 1) for the gene.
#'
#' @return A scatter plot representing for each TCGA sample, gene expression 
#'     and mean methylation values of probe(s) located in its promoter region 
#'     (defined as 1000 nucleotides upstream TSS and 200 nucleotides downstream
#'     TSS by default). If return = TRUE, methylation and expression values 
#'     are returned in a tibble instead. 
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_point ggtitle xlim scale_colour_manual
#' @importFrom stats cor.test quantile
#' @importFrom rlang .data
#'
#' @examples
#' TCGA_methylation_expression_correlation("LUAD", gene = "TDRD1")
#' TCGA_methylation_expression_correlation("all", gene = "MAGEA1", 
#'  min_probe_number = 1, return = FALSE)
TCGA_methylation_expression_correlation <- function(
    tumor = "all",
    gene = NULL,
    nt_up = 1000,
    nt_down = 200,
    min_probe_number = 3,
    include_normal_tissues = FALSE,
    return = FALSE) {

    methylation_expression <- prepare_TCGA_methylation_expression(
      tumor = tumor, 
      gene = gene,
      nt_up = nt_up,
      nt_down = nt_down,
      include_normal_tissues = include_normal_tissues)
    
    methylation_expression <- methylation_expression[
      methylation_expression$probe_number >= min_probe_number, ]
    methylation_expression <- methylation_expression[
      !is.na(methylation_expression$met), ]
    
    if (nrow(methylation_expression) == 0) stop("less than ",  
                                                min_probe_number, 
                                                " probes in selected region")
    
    ## Evaluate correlation only if the gene is expressed (TPM >= 1) 
    ## in at least 1% of the samples
    if (quantile(methylation_expression$TPM, 0.99) < 1) {
      message("Too few positive samples to estimate a correlation for ", gene)
      cor <- NA
    }
    
    if (!is.na(cor)) cor <- cor.test(methylation_expression$met,
                                     log1p(methylation_expression$TPM))$estimate
  
    p <- ggplot(
            as_tibble(methylation_expression[order(methylation_expression$type,
                decreasing = TRUE), ]),
            aes(x = .data$met, y = log1p(.data$TPM))) +
            geom_point(alpha = 0.8, aes(color = .data$type,
                                        shape = .data$tissue)) +
            ggtitle(paste0(gene, "(Pearson's corr = ", round(cor, 2), ")")) +
            scale_colour_manual(values = TCGA_colors) +
            xlim(0, 1)

    if (return) {
        return(methylation_expression)
    }

    return(p)
}

