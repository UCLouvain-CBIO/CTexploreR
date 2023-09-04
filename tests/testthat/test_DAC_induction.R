test_that("DAC_induction() works", {
  
  ## Check returns a matrix of double
  res <- DAC_induction(c("MAGEA1", "MAGEA3"), return = TRUE)
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  
  ## n valid genes in input returns a matrix of n expected rownames
  expect_equal(nrow(res), 2) 
  expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
  
  res <- DAC_induction("MAGEA1", return = TRUE)
  expect_equal(nrow(res), 1) 
  
  ## Check fake name gives nothing and warning 
  expect_equal(nrow(DAC_induction("", return = TRUE)), 0) 
  expect_warning(DAC_induction(""), "names invalid")
  
  ## Test the multimapping argument
  res_mm <- DAC_induction("MAGEA2", return = TRUE)
  res_no_mm <- DAC_induction("MAGEA2", multimapping = FALSE, return = TRUE)
  expect_gt(sum(res_mm), sum(res_no_mm)) 
  
  ## Test that the function returns a heatmap
  res <- DAC_induction(c("MAGEA1", "MAGEA3"))
  expect_s4_class(res, "Heatmap") 
  
  ## Check no variable gives CT_genes
  res <- DAC_induction(return = TRUE)
  expect_identical(rownames(res),
                   CT_genes$external_gene_name)
  expect_equal(nrow(res), nrow(CT_genes))
  
   
})
