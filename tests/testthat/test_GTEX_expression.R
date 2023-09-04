test_that("GTEX_expression() works", {
    
    ## returns a matrix of double
    res <- CTexploreR::GTEX_expression(c("MAGEA1", "MAGEA3"), return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## n valid genes in input returns a matrix of n expected rownames
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
        
    res <- CTexploreR::GTEX_expression("MAGEA1", return = TRUE)
    expect_equal(nrow(res), 1) 
    
    res <- CTexploreR::GTEX_expression("", return = TRUE)
    expect_equal(nrow(res), 0) 
    expect_warning(CTexploreR::GTEX_expression(my_genes), "names invalid")
    
    ## Test the "log_TPM" units argument
    res_in_TPM <- CTexploreR::GTEX_expression("MAGEA1", return = TRUE)
    res_in_log <- CTexploreR::GTEX_expression("MAGEA1", 
                                               units = "log_TPM", return = TRUE)
    expect_equal(res_in_log[, "Testis"], log1p(res_in_TPM[, "Testis"])) 
    
    ## Test that the function returns a heatmap
    res <- CTexploreR::GTEX_expression(c("MAGEA1", "MAGEA3"))
    expect_s4_class(res, "Heatmap")
})
