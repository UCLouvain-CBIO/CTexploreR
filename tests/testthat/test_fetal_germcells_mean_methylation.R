test_that("fetal_germcells_mean_methylation() works", {
  
  ## returns a SummarizedExperiment when values_only is set to TRUE
  res <- fetal_germcells_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4"), 
    values_only = TRUE)
  expect_s4_class(res, "SummarizedExperiment")  
  expect_equal(nrow(res), 3)

  
  ## no valid gene in input returns an empty matrix
  expect_warning(res_no_gene <- fetal_germcells_mean_methylation(
    genes = "xxx", values_only = TRUE), "valid types")
  expect_equal(nrow(res_no_gene), 0)
  
  ## Test that the function returns a heatmap by default
  ## res <- fetal_germcells_mean_methylation(
  ##  genes = c("MAGEA1", "MAGEA2", "MAGEA4", "TDRD1"))
  ## expect_s4_class(res, "Heatmap")
})
