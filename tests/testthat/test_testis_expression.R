test_that("testis_expression() works", {
    
    ## returns a SingleCellExperiment
    res <- testis_expression(cells = "germ_cells", 
                             genes = c("MAGEA1", "MAGEA3"), 
                             values_only = TRUE)
    expect_s4_class(res, "SingleCellExperiment")
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
    
    ## selects the expected cell types
    ## works with only one gene in input
    my_cells <- c("Spermatogonia", "Sertoli")
    res <- testis_expression(cells = my_cells, genes = "MAGEA1", 
                             values_only = TRUE)
    expect_equal(nrow(res), 1) 
    expect_true(all(unique(res$type) %in% my_cells))
    
    ## selects correctly germ_cells
    my_cells <- "germ_cells"
    res <- testis_expression(cells = my_cells, genes = "MAGEA1", 
                             values_only = TRUE)
    germ_cells <- c("SSC", "Spermatogonia", "Early_spermatocyte",
                    "Late_spermatocyte", "Round_spermatid", 
                    "Elongated_spermatid", "Sperm1", "Sperm2")
    expect_true(all(unique(res$type) %in% germ_cells))
    
    ## selects correctly somatic_cells
    my_cells <- "somatic_cells"
    res <- testis_expression(cells = my_cells, genes = "MAGEA1", 
                             values_only = TRUE)
    somatic_cells <- c("Macrophage", "Endothelial",
                       "Myoid", "Sertoli", "Leydig")
    expect_true(all(unique(res$type) %in% somatic_cells))
    
    ## Test that the function returns a heatmap
    ## res <- testis_expression(cells = germ_cells, 
    ##                         genes = c("MAGEA1", "MAGEA3", "MAGEA4", "MAGEC1"))
    ## expect_s4_class(res, "Heatmap")
    ## vdiffr::expect_doppelganger("testis_expression_on_MAGE", fig = res)
})
