test_that("testis_expression() works", {
    
    ## returns a SingleCellExperiment
    res <- testis_expression(cells = "germ_cells", 
                             genes = c("MAGEA1", "MAGEA3"), 
                             return = TRUE)
    expect_s4_class(res, "SingleCellExperiment")
    
    ## n valid genes in input returns a matrix of n expected rownames
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
    
    res <- testis_expression(cells = "germ_cells", 
                             genes = "MAGEA4", 
                             return = TRUE)
    expect_equal(nrow(res), 1) 

    res <- testis_expression(cells = "germ_cells", 
                             genes = "xxx", 
                             return = TRUE)
    expect_equal(nrow(res), 0)
    expect_warning(testis_expression(cells = "germ_cells", genes = "xxx", 
                                     return = TRUE), "names invalid")
    
    ## selects the expected cell types
    my_cells <- c("Spermatogonia", "Sertoli")
    res <- testis_expression(cells = my_cells, genes = "MAGEA1", return = TRUE)
    x <- colData(CTdata::testis_sce())
    exp_cells <- rownames(x[x$type %in% my_cells, ])
    expect_equal(sort(colnames(res)), sort(exp_cells))
    
    ## selects correctly germ_cells
    my_cells <- "germ_cells"
    res <- testis_expression(cells = my_cells, genes = "MAGEA1", return = TRUE)
    germ_cells <- c("SSC", "Spermatogonia", "Early_spermatocyte",
                    "Late_spermatocyte", "Round_spermatid", 
                    "Elongated_spermatid", "Sperm1", "Sperm2")
    exp_cells <- rownames(x[x$type %in% germ_cells, ])
    expect_equal(sort(colnames(res)), sort(exp_cells))
    
    ## selects correctly somatic_cells
    my_cells <- "somatic_cells"
    res <- testis_expression(cells = my_cells, genes = "MAGEA1", return = TRUE)
    somatic_cells <- c("Macrophage", "Endothelial",
                       "Myoid", "Sertoli", "Leydig")
    exp_cells <- rownames(x[x$type %in% somatic_cells, ])
    expect_equal(sort(colnames(res)), sort(exp_cells))
    
    ## Test that the function returns a heatmap
    res <- testis_expression(cells = my_cells, genes = "MAGEA1")
    expect_s4_class(res, "Heatmap")
})
