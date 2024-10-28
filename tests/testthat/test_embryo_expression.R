test_that("embryo_expression() works", {
  
  ## returns a SingleCellExperiment
  ## n valid genes in input returns a SingleCellExperiment of n expected rownames
  ## returns a warning when an invalid gene is entered
  ## works with only one gene in input
  ## error if no dataset is specified
  expect_error(embryo_expression("MAGEA1"))
  expect_warning(res <- embryo_expression(c("MAGEA1", "xxx"),
                                          dataset = "Petropoulos",
                                          values_only = TRUE), "names invalid")
  expect_s4_class(res, "SingleCellExperiment")  
  expect_equal(nrow(res), 1) 
  expect_identical(rownames(res), "MAGEA1")
  

  ## Test that the function returns a heatmap by default
  ## Test that returns the expected heatmap
  ## res <- embryo_expression(c("MAGEA1", "MAGEA3", "MAGEA4"),
  ##                         dataset = "Petropoulos")
  ## expect_s4_class(res, "Heatmap")
  
  ## Same with Zhu
  expect_warning(res <- embryo_expression(c("MAGEA1", "xxx"),
                                          dataset = "Zhu",
                                          values_only = TRUE), "names invalid")
  expect_s4_class(res, "SingleCellExperiment")  
  expect_equal(nrow(res), 1) 
  expect_identical(rownames(res), "MAGEA1")
  
})
  