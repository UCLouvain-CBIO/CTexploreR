test_that("CT_correlated_genes() works", {

    ## returns a tibble when Return is set to TRUE
    ## Correlations corresponding to all genes vs specified CT gene
    ## CT_genes column specifies if genes are CT genes
    res <- CT_correlated_genes(c("MAGEA3"), values_only = TRUE)
    expect_true(inherits(res, "data.frame"))
    expect_equal(nrow(res), 24473)
    expect_identical(res$external_gene_name[1], "MAGEA3")
    expect_true(all(res$external_gene_name[res$CT_gene] %in%
                      CT_genes$external_gene_name))
    expect_true(!all(res$external_gene_name[!res$CT_gene] %in%
                       CT_genes$external_gene_name))

    ## returns an error if no specified gene or if the gene is not a CT gene
    expect_error(CT_correlated_genes(gene = "BRCA1"), "CT gene")
})
