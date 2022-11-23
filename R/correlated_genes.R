#' Gene correlations in CCLE cancer cell lines
#'
#' @description A function that uses expression data from CCLE cell lines and
#' highlights genes correlated (or anti-correlated) with specified CT gene.
#'
#' @param corr_matrix CCLE_correlation_matrix
#'
#' @param gene CT gene selected
#'
#' @param corr_thr Genes with an absolute correlation coefficient higher
#' than this threshold will be highlighted (default = 0.5)
#'
#' @return A plot where each dots represent the correlation coefficients between
#' genes and the specified CT gene (entered as input). All correlations coefficients
#' are invisibly returned.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom ggrepel geom_text_repel
#'
#' @examples
#' correlated_genes(CCLE_correlation_matrix, gene = "MAGEA3")
#' correlated_genes(CCLE_correlation_matrix, "TDRD1", 0.3)
correlated_genes <- function(corr_matrix, gene, corr_thr = 0.5){

  if (missing(corr_matrix)) {
    stop("Correlation matrix be specified!")
  }

  if (missing(gene)) {
    stop("Gene name be specified!")
  }

  if (!gene %in% CT_genes$external_gene_name) {
    stop("Gene must be a CT gene!")
  }

  tested_ref <- rowData(CCLE_data)[rowData(CCLE_data)$external_gene_name == gene, "ensembl_gene_id"]

  tmp <- data.frame(ensembl_gene_id = names(corr_matrix[tested_ref, ]),
                    corr = corr_matrix[tested_ref, ],
                    external_gene_name =rowData(CCLE_data)[names(corr_matrix[tested_ref, ]), "external_gene_name"])
  tmp$CT_gene <- ifelse(tmp$external_gene_name %in% CT_genes$external_gene_name, "CTgene", "not_CTgene")
  tmp <- tmp[order(tmp$corr, decreasing = TRUE), ]
  highly_correlated <- tmp[(!is.na(tmp$corr) & (tmp$corr > corr_thr | tmp$corr < -corr_thr)), "external_gene_name"]

  p <- ggplot(tmp[tmp$external_gene_name %in% highly_correlated, ],
              aes(x = gene, y = corr, label = external_gene_name)) +
    geom_jitter(alpha = 0.5, color = "blue", position = position_jitter(height = 0, seed = 1)) +
    geom_text_repel(position = position_jitter(height = 0, seed = 1), size = 2.5, max.overlaps = 20) +
    geom_jitter(data = tmp[!tmp$external_gene_name %in% highly_correlated, ],
                aes(x = gene, y = corr), alpha = 0.3) +
    ggtitle(paste0("Genes correlated with ", gene)) +
    geom_hline(yintercept = -corr_thr, linetype = "dashed", size = 0.5, color = "blue") +
    geom_hline(yintercept = 0, size = 0.5, color = "blue") +
    geom_hline(yintercept = corr_thr, linetype = "dashed", size = 0.5,  color = "blue") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("Correlation coefficient") +
    ylim(-1, 1)

  suppressWarnings(print(p))
  invisible(tmp)
}
