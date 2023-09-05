test_that("HPA_cell_type_expression() works", {
    
    ## returns a SingleCellExperiment when Return is set to TRUE
    res <- HPA_cell_type_expression("MAGEA1", return = TRUE)
    expect_s4_class(res, "SingleCellExperiment")
    
    ## Test that the function returns a heatmap by default
    res <- HPA_cell_type_expression("MAGEA1")
    expect_s4_class(res, "Heatmap")
    
    ## n valid genes in input returns a matrix of n expected rownames
    res <- HPA_cell_type_expression(c("MAGEA1", "MAGEA3"), return = TRUE)
    expect_equal(nrow(res), 2) 
    returned_genes <- pull(CT_genes[CT_genes$ensembl_gene_id %in% rownames(res), 
                               "external_gene_name"])
    expect_identical(sort(returned_genes), sort(c("MAGEA1", "MAGEA3")))
    
    ## works with only one gene in input
    res <- HPA_cell_type_expression("MAGEA1", return = TRUE)
    expect_equal(nrow(res), 1)
    
    ## works but returns a warning when no valid gene is entered
    res <- HPA_cell_type_expression("", return = TRUE)
    expect_equal(nrow(res), 0)
    expect_warning(HPA_cell_type_expression("", return = TRUE), "names invalid")
})
