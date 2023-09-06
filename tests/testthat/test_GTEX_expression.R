test_that("GTEX_expression() works", {
    
    ## returns a matrix of double
    res <- GTEX_expression(c("MAGEA1", "MAGEA3"), return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## n valid genes in input returns a matrix of n expected rownames
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
        
    ## works with only one gene in input
    res <- GTEX_expression("MAGEA1", return = TRUE)
    expect_equal(nrow(res), 1) 
    
    ## works but returns a warning when no valid gene is entered
    res <- GTEX_expression("", return = TRUE)
    expect_equal(nrow(res), 0) 
    expect_warning(GTEX_expression(""), "names invalid")
    
    ## Test the "log_TPM" units argument
    res_in_TPM <- GTEX_expression("MAGEA1", return = TRUE)
    res_in_log <- GTEX_expression("MAGEA1", units = "log_TPM", return = TRUE)
    expect_equal(res_in_log[, "Testis"], log1p(res_in_TPM[, "Testis"])) 
    
    ## Test that the function returns a heatmap
    res <- GTEX_expression(c("MAGEA1", "MAGEA3"))
    expect_s4_class(res, "Heatmap")
    
    ## Test that the function returns the expected heatmap
    res <- GTEX_expression(c("MAGEA1", "MAGEA3", "MAGEA4"))
    vdiffr::expect_doppelganger("GTEX_expression on MAGE", fig = res)
})
