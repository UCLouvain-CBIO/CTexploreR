test_that("GTEX_expression() works", {
  
  ## returns a matrix of double
  ## n valid genes in input returns a matrix of n expected rownames
  ## returns a warning when an invalid gene is entered
  ## Test that the function returns the expected matrix
  ## works with only one gene in input
  expect_warning(res <- GTEX_expression(c("MAGEA1", "xxx"), 
                                        return = TRUE), "names invalid")
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 1) 
  expect_identical(rownames(res), "MAGEA1")
  #saveRDS(res, test_path("fixtures", "GTEX_expression_on_MAGE.rds"))
  exp_res <- readRDS(test_path("fixtures", "GTEX_expression_on_MAGE.rds"))
  expect_equal(exp_res, res)
  
  ## Test the "log_TPM" units argument
  res_in_log <- GTEX_expression(c("MAGEA1"), 
                                units = "log_TPM", 
                                return = TRUE)
  expect_equal(res_in_log[, "Testis"], log1p(res[, "Testis"])) 
  
  ## Test that the function returns a heatmap by default
  ## Test that returns the expected heatmap
  ## res <- GTEX_expression(c("MAGEA1", "MAGEA3", "MAGEA4"))
  ## expect_s4_class(res, "Heatmap")
  ## vdiffr::expect_doppelganger("GTEX_expression_on_MAGE", fig = res)
    
})

