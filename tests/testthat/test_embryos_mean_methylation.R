test_that("embryos_mean_methylation() works", {
  
  ## returns a RangedSummarizedExperiment when values_only is set to TRUE
  res <- embryos_mean_methylation(
    genes = c("MAGEA1", "MAGEA2", "MAGEA4"), 
    values_only = TRUE)
  expect_s4_class(res, "RangedSummarizedExperiment")  
  expect_equal(nrow(res), 3)
  
  ## selects the expected cell types
  my_stages <- c("Morula", "Blastocyst")
  res <- embryos_mean_methylation(stage = my_stages, genes = "MAGEA1",
                                  values_only = TRUE)
  expect_equal(nrow(res), 1) 
  expect_true(all(unique(colData(res)$Stage) %in% my_stages))
  
  ## no valid gene in input returns an empty matrix
  expect_warning(res_no_gene <- embryos_mean_methylation(
    genes = "xxx", values_only = TRUE), "valid types")
  expect_equal(nrow(res_no_gene), 0)
  
  ## Test that the function returns a heatmap by default
  ## res <- embryos_mean_methylation(
  ##  genes = c("MAGEA1", "MAGEA2", "MAGEA4", "TDRD1"))
  ## expect_s4_class(res, "Heatmap")
})
