test_that("DAC_induction() works", {
  
  ## Check returns a matrix of double
  ## n valid genes in input returns a matrix of n expected rownames
  expect_warning(res <- DAC_induction(c("MAGEA1", "MAGEA3", "xxx"), 
                                      values_only = TRUE), "names invalid")
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 2) 
  expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
    
  ## Check fake name gives nothing and warning 
  expect_warning(expect_equal(nrow(DAC_induction("", values_only = TRUE)), 0), 
                 "names invalid")
  
  ## Test the multimapping argument
  ## Works wit one gene only
  res_mm <- DAC_induction("MAGEA3", values_only = TRUE)
  res_no_mm <- DAC_induction("MAGEA3", multimapping = FALSE, values_only = TRUE)
  expect_gt(sum(res_mm), sum(res_no_mm)) 
  
  ## Check that only strict CT_genes are returned if no gene is specified
  res <- DAC_induction(values_only = TRUE)
  expect_equal(nrow(res), 146)
  
  ## Test that the function returns a heatmap
  ## res <- DAC_induction(c("MAGEA1", "MAGEA3"))
  ## expect_s4_class(res, "Heatmap") 
  ## vdiffr::expect_doppelganger("DAC_induction_MAGEA1_MAGEA3", fig = res)
})
