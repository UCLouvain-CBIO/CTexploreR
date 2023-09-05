test_that("CCLE_expression() works", {
    
    ## returns a matrix of double when Return is set to TRUE
    res <- CCLE_expression("MAGEA1", type = "Skin", return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## Test that the function returns a heatmap by default
    res <- CCLE_expression("MAGEA1", "Skin")
    expect_s4_class(res, "Heatmap")
    
    ## works with only one gene in input
    res <- CCLE_expression("MAGEA1", type = "Skin", return = TRUE)
    expect_equal(nrow(res), 1) 
    
    ## n valid genes in input returns a matrix of n expected rownames
    res <- CCLE_expression(c("MAGEA1", "MAGEA3"), type = "Skin", return = TRUE)
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
        
    ## works but returns a warning when no valid gene is entered
    expect_equal(nrow(CCLE_expression("", type = "Skin", return = TRUE)), 0)
    expect_warning(CCLE_expression("", type = "Skin"), "names invalid")
    
    ## No valid tumor type returns an error
    expect_error(CCLE_expression(genes = "MAGEA1"), "No valid")
    
    ## selects the expected cell lines
    my_types <- c("Skin", "Colorectal")
    res <- CCLE_expression(genes = "MAGEA1", type = my_types, return = TRUE)
    x <- colData(CCLE_data())
    exp_cells <- rownames(x[x$type %in% my_types, ])
    expect_equal(sort(colnames(res)), sort(exp_cells))
    
    ## tumor type is case insensitive
    my_types <- c("sKIN")
    res <- CCLE_expression(genes = "MAGEA1", type = my_types, return = TRUE)
    expect_equal(table(x$type == "Skin")[["TRUE"]], ncol(res))
    
    ## Test the "log_TPM" units argument
    res_in_TPM <- CCLE_expression("MAGEA1", "Skin", return = TRUE)
    res_in_log <- CCLE_expression("MAGEA1", "Skin", units = "log_TPM", 
                                  return = TRUE)
    expect_equal(res_in_log[, 1], log1p(res_in_TPM[, 1])) 
    
})
