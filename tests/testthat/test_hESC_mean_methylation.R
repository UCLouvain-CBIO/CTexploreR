test_that("hESC_mean_methylation() works", {
  
  ## returns a SummarizedExperiment when values_only is set to TRUE
  ## the na.omit parameter set to FALSE keeps genes with missing values
  res <- hESC_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4"), 
    na.omit = FALSE, values_only = TRUE)
  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(nrow(res), 3)
  
  ## the na.omit parameter set to TRUE removes genes with missing values
  res_na_omit <- hESC_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4"), 
    na.omit = TRUE, values_only = TRUE)
  expect_identical(na.omit(res), res_na_omit)
  
  ## no valid gene in input returns an empty matrix
  expect_warning(res_no_gene <- hESC_mean_methylation(
    genes = "xxx", values_only = TRUE), "valid types")
  expect_equal(nrow(res_no_gene), 0)
  
  ## Test that the function returns a heatmap by default
  ## res <- hESC_mean_methylation(
  ##  genes = c("MAGEA1", "MAGEA2", "MAGEA4", "TDRD1"))
  ## expect_s4_class(res, "Heatmap")

})