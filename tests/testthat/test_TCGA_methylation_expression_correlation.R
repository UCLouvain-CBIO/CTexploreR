test_that("TCGA_methylation_expression_correlation() works", {
    
    ## returns a tibble when Return is set to TRUE
    res <- TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                   gene = "TDRD1",
                                                   return = TRUE)
    expect_true(inherits(res, "data.frame"))
    
    ## same output as prepare_TCGA_methylation_expression if Return is TRUE
    expected_res <- prepare_TCGA_methylation_expression(tumor = "LUAD", 
                                                       gene = "TDRD1")
    expect_identical(res, expected_res)
    
    ## returns a plot when Return is set to TRUE
    res <- TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                   gene = "TDRD1",
                                                   return = FALSE)
    expect_s3_class(res, "ggplot")
    
    ## results are filtered according to the minimum number of probe to include
    res <- TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                   gene = "TDRD1",
                                                   min_probe_number = 6,
                                                   return = TRUE)
    expect_true(all(res$probe_number >= 6)) 
    
    
    ## returns an error if not enough probes to evaluate the mean methylation
    expect_error(TCGA_methylation_expression_correlation(tumor = "LUAD", 
                                                         gene = "MAGEA1",
                                                         min_probe_number = 3,
                                                         return = TRUE), 
                 "probes")

})
