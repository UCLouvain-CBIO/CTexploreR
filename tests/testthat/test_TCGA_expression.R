test_that("TCGA_expression() works", {
    
    ## returns a matrix of double
    res <- TCGA_expression("LUAD", "MAGEA1", return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## n valid genes in input returns a matrix of n expected rownames
    expect_equal(nrow(res), 1) 
    
    res <- TCGA_expression("LUAD", c("MAGEA1", "MAGEA3"), return = TRUE)
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
        
    res <- TCGA_expression("LUAD", "", return = TRUE)
    expect_equal(nrow(res), 0)
    expect_warning(TCGA_expression("LUAD", "", return = TRUE), "names invalid")
    
    ## selects the expected tumors
    tumor_types <- c("SKCM")
    res <- TCGA_expression(tumor = tumor_types, genes = "MAGEA1", return = TRUE)
    x <- colData(TCGA_TPM())
    exp_samples <- rownames(x[x$project_id %in% paste0("TCGA-", tumor_types), ])
    expect_equal(sort(colnames(res)), sort(exp_samples))
    
    ## Peritumoral samples are displayed only when a single type of tumor asked
    tumor_types <- c("SKCM", "LUSC")
    res <- TCGA_expression(tumor = tumor_types, genes = "MAGEA1", return = TRUE)
    exp_samples <- rownames(x[x$project_id %in% paste0("TCGA-", tumor_types) &
                                x$shortLetterCode != "NT", ])
    expect_equal(sort(colnames(res)), sort(exp_samples))
    
    ## No valid tumor type returns an error
    expect_error(TCGA_expression("", "MAGEA1"), "No valid")
    
    ## Test the "log_TPM" units argument
    res_in_TPM <- TCGA_expression("LUAD", "MAGEA1", return = TRUE)
    res_in_log <- TCGA_expression("LUAD", "MAGEA1", units = "log_TPM", 
                                  return = TRUE)
    expect_equal(res_in_log[, 1], log1p(res_in_TPM[, 1])) 
    
    ## Test that the function returns a heatmap
    res <- TCGA_expression("LUAD", "MAGEA1")
    expect_s4_class(res, "Heatmap")
})



  