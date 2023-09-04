test_that("CCLE_expression() works", {
    
    ## returns a matrix of double
    res <- CTexploreR:::CCLE_expression("MAGEA1", type = "Skin", return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## n valid genes in input returns a matrix of n expected rownames
    expect_equal(nrow(res), 1) 
    
    res <- CTexploreR:::CCLE_expression(c("MAGEA1", "MAGEA3"), 
                                        type = "Skin", return = TRUE)
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
        
    res <- CTexploreR:::CCLE_expression("", type = "Skin", return = TRUE)
    expect_equal(nrow(res), 0)
    expect_warning(CTexploreR:::CCLE_expression(my_genes, type = "Skin"), 
                   "names invalid")
    
    ## selects the expected cell lines
    my_types <- c("Skin", "Colorectal")
    res <- CTexploreR:::CCLE_expression(genes = "MAGEA1", 
                                        type = my_types, 
                                        return = TRUE)
    x <- SummarizedExperiment:::colData(CTdata:::CCLE_data())
    exp_cells <- rownames(x[x$type %in% my_types, ])
    expect_equal(sort(colnames(res)), sort(exp_cells))
    
    ## tumor type is case insensitive
    my_types <- c("sKIN")
    res <- CTexploreR:::CCLE_expression(genes = "MAGEA1", 
                                        type = my_types, return = TRUE)
    expect_equal(table(x$type == "Skin")[["TRUE"]], ncol(res))
    
    ## No valid tumor type returns an error
    expect_error(CTexploreR:::CCLE_expression(genes = "MAGEA1"), "No valid")

    ## Test the "log_TPM" units argument
    res_in_TPM <- CTexploreR:::CCLE_expression("MAGEA1", "Skin", return = TRUE)
    res_in_log <- CTexploreR:::CCLE_expression("MAGEA1", "Skin", 
                                               units = "log_TPM", return = TRUE)
    expect_equal(res_in_log[, 1], log1p(res_in_TPM[, 1])) 
    
    ## Test that the function returns a heatmap
    res <- CTexploreR:::CCLE_expression("MAGEA1", "Skin")
    expect_s4_class(res, "Heatmap")
})
