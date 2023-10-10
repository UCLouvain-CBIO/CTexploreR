test_that("normal_tissues_mean_methylation() works", {
  
  ## returns a matrix of double when Return is set to TRUE
  ## the na.omit parameter set to FALSE keeps genes with missing values
  res <- normal_tissues_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4"), 
    na.omit = FALSE, return = TRUE)
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 3)
  
  ## the na.omit parameter set to TRUE removes genes with missing values
  res_na_omit <- normal_tissues_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4"), 
    na.omit = TRUE, return = TRUE)
  expect_identical(na.omit(res), res_na_omit)
  
  ## no valid gene in input returns an empty matrix
  expect_warning(res_no_gene <- normal_tissues_mean_methylation(
    genes = "xxx", return = TRUE), "valid types")
  expect_equal(nrow(res_no_gene), 0)

  ## Test that the function returns a heatmap by default
  res <- normal_tissues_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4", "TDRD1"))
  expect_s4_class(res, "Heatmap")
  vdiffr::expect_doppelganger("normal_tissues_mean_methylation_on_few_genes", fig = res)

  
})
