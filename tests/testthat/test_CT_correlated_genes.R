test_that("CT_correlated_genes() works", {
    
    ## returns a tibble when Return is set to TRUE
    res <- CT_correlated_genes(gene = "MAGEA3", return = TRUE)
    expect_true(inherits(res, "data.frame"))
    
    ## returns a plot when Return is set to TRUE
    res <- CT_correlated_genes(gene = "MAGEA3", return = FALSE)
    expect_s3_class(res, "ggplot")
    
    ## CT_genes column specifies if genes are CT genes
    res <- CT_correlated_genes(gene = "MAGEA3", return = TRUE)
    expect_true(all(res$external_gene_name[res$CT_gene] %in% 
                      CT_genes$external_gene_name))
    expect_true(!all(res$external_gene_name[!res$CT_gene] %in% 
                      CT_genes$external_gene_name))
    
    ## returns an error if no specified gene or if the gene is not a CT gene
    expect_error(CT_correlated_genes(), "Gene")
    expect_error(CT_correlated_genes(gene = "BRCA1"), "CT gene")
    
})
