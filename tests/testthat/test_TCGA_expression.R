test_that("TCGA_expression() works", {
  
  ## returns a matrix of double
  ## returns a warning when an invalid gene is entered
  ## n valid genes in input returns a matrix of n expected rownames
  expect_warning(res <- TCGA_expression("LUAD", c("MAGEA1", "MAGEA3", "xxx"), 
                                        return = TRUE), "names invalid")
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 2) 
  expect_identical(sort(rownames(res)), sort(c("MAGEA1", "MAGEA3")))
  
  ## Test the "log_TPM" units argument
  ## returns a warning when an invalid tumor type is entered
  expect_warning(res_in_log <- TCGA_expression(c("LUAD", "lung"), 
                                               c("MAGEA1", "MAGEA3"), 
                                units = "log_TPM", return = TRUE), 
                 "names invalid")
  expect_equal(res_in_log[, 1:3], log1p(res[, 1:3])) 

  ## selects the expected number of LUAD samples
  ## x <- colData(TCGA_TPM())
  ## exp_samples <- rownames(x[x$project_id %in% "TCGA-LUAD", ])
  ## exp_number_of_LUAD_samples <- length(exp_samples)
  exp_number_of_LUAD_samples <- 598
  expect_equal(dim(res)[2], exp_number_of_LUAD_samples)
  
  ## Peritumoral samples are returned when a single type of tumor asked
  expect_true(length(grep("TCGA-\\d{2}-\\d{4}-11", x = colnames(res))) > 0)
  
  ## Test that the function returns a heatmap by default
  ## res <- TCGA_expression(c("LUAD", "SKCM"), "MAGEA1")
  ## expect_s4_class(res, "Heatmap")
  ## vdiffr::expect_doppelganger("TCGA_expression_on_MAGEA1", 
  ##                             fig = res)
})



