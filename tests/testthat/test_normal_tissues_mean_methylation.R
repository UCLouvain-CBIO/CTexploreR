test_that("normal_tissues_mean_methylation() works", {
    
    ## returns a matrix of double
    res <- normal_tissues_mean_methylation(
      genes = c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"), 
      na.omit = FALSE, return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
  
    ## n valid genes in input returns a matrix of n expected rownames
    expect_equal(nrow(res), 4) 
    
    ## no valid gene in input returns an empty matrix
    res_no_gene <- normal_tissues_mean_methylation(
      genes = "xxx", 
      na.omit = FALSE, return = TRUE)
    expect_equal(nrow(res_no_gene), 0)
    
    ## the na.omit parameter set to TRUE removes gens with missing values
    res_na_omit <- normal_tissues_mean_methylation(
      genes = c("MAGEA1", "MAGEA2", "MAGEA3", "MAGEA4"), 
      na.omit = TRUE, return = TRUE)
    expect_identical(na.omit(res), res_na_omit)

    ## collects the correct data
    res <- normal_tissues_mean_methylation(genes = "MAGEA1", return = TRUE)
    expect_identical(res, assay(CT_mean_methylation_in_tissues()["MAGEA1",]))
    
    ## Test that the function returns a heatmap
    res <- normal_tissues_mean_methylation()
    expect_s4_class(res, "Heatmap")
})
