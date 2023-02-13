#' Gene correlations in CCLE cancer cell lines
#'
#' @description A function that uses expression data from CCLE cell lines and
#' highlights genes correlated (or anti-correlated) with specified CT gene.
#' Genes with a correlation coefficient above threshold are colored in red if
#' they are CT genes or in blue, if not.
#'
#' @param gene CT gene selected
#'
#' @param corr_thr Genes with an absolute correlation coefficient (Pearson)
#' higher than this threshold will be highlighted (default = 0.5)
#'
#' @param return Boolean (FALSE by default). If set to TRUE, the function will
#' return the correlation coefficients with all genes instead of the plot.
#'
#' @return A plot where each dots represent the correlation coefficients (Pearson)
#' between genes and the specified CT gene (entered as input). Genes with a
#' correlation coefficient above threshold are colored in red if they are CT
#' genes or in blue, if not. If return = TRUE, all correlations coefficients are
#' returned instead.
#'
#' @export
#'
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 geom_jitter position_jitter scale_colour_manual geom_hline theme element_blank ylab ylim
#' @importFrom rlang .data
#'
#' @examples
#' correlated_genes(gene = "MAGEA3")
#' correlated_genes("TDRD1", 0.3)
correlated_genes <- function(gene, corr_thr = 0.5,
                             return = FALSE) {

    corr_matrix <- CTdata::CCLE_correlation_matrix()
    CT_genes <- CTdata::CT_genes()

    if (missing(gene)) {
        stop("Gene name be specified!")
    }

    if (!gene %in% CT_genes$external_gene_name) {
        stop("Gene must be a CT gene!")
    }

    tested_ref <- rownames(CCLE_data[rowData(CCLE_data)$external_gene_name == gene, ])

    tmp <- data.frame(corr = corr_matrix[tested_ref, ],
                      external_gene_name = rowData(CCLE_data)[names(corr_matrix[tested_ref, ]),
                                                              "external_gene_name"])
    tmp$CT_gene <- ifelse(tmp$external_gene_name %in% CT_genes$external_gene_name,
                          TRUE, FALSE)
    tmp <- tmp[order(tmp$corr, decreasing = TRUE), ]
    highly_correlated <-
        tmp[(!is.na(tmp$corr) & (tmp$corr > corr_thr | tmp$corr < -corr_thr)),
            "external_gene_name"]

    p <- ggplot(tmp[tmp$external_gene_name %in% highly_correlated, ],
                aes(x = gene, y = .data$corr, 
                    label = .data$external_gene_name)) +
        geom_jitter(aes(color = .data$CT_gene), alpha = 0.5,
                    position = position_jitter(height = 0, seed = 1)) +
        scale_colour_manual(limits = c(TRUE, FALSE),
                            values = c("red", "deepskyblue")) +
        geom_text_repel(position = position_jitter(height = 0, seed = 1),
                        size = 2.5, max.overlaps = 20) +
        geom_jitter(data = tmp[!tmp$external_gene_name %in% highly_correlated, ],
                    aes(x = gene, y = .data$corr), alpha = 0.3) +
        ggtitle(paste0("Genes correlated with ", gene)) +
        geom_hline(yintercept = -corr_thr, linetype = "dashed",
                   linewidth = 0.5, color = "blue") +
        geom_hline(yintercept = 0, linewidth = 0.5, color = "blue") +
        geom_hline(yintercept = corr_thr, linetype = "dashed",
                   linewidth = 0.5,  color = "blue") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom") +
        ylab("Correlation coefficient") +
        ylim(-1, 1)


    if (return)
        return(tmp)

    suppressWarnings(print(p))
}
