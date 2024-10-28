test_that("oocytes_expression() works", {
  
  ## returns a SingleCellExperiment
  ## n valid genes in input returns a SingleCellExperiment of n expected rownames
  ## returns a warning when an invalid gene is entered
  ## works with only one gene in input
  
  expect_warning(res <- oocytes_expression(c("MAGEA1", "xxx"),
                                                   values_only = TRUE), "names invalid")
  expect_s4_class(res, "SingleCellExperiment")  
  expect_equal(nrow(res), 1) 
  expect_identical(rownames(res), "MAGEA1")
  
  
  ## Test that the function returns a heatmap by default
  ## Test that returns the expected heatmap
  ## res <- oocytes_expression(c("MAGEA1", "MAGEA3", "MAGEA4"))
  ## expect_s4_class(res, "Heatmap")
  
})
