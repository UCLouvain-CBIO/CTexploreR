test_that("prepare_TCGA_methylation_expression_() works", {
  
  ## returns a tibble
  res <- prepare_TCGA_methylation_expression(tumor = "all", 
                                             gene = "TDRD1")
  expect_true(inherits(res, "data.frame"))

  ## returns an error if no valid tumor type is entered. 
  expect_error(prepare_TCGA_methylation_expression(
    tumor = "xxx", 
    gene = "TDRD1"), 
    "No valid")
  
  ## returns values for all tumors when tumor parameter is set to "all"
  expect_true(all(unique(res$type) %in% c("SKCM", "LUAD", "LUSC", "COAD", 
                                          "ESCA", "BRCA", "HNSC")))
  
  ## returns values for the specified tumor type 
  res <- prepare_TCGA_methylation_expression(tumor = "LUSC", 
                                             gene = "TDRD1")
  expect_true(all(unique(res$type) == "LUSC"))
  
  ## returns an error if no valid gene name entered  
  expect_error(prepare_TCGA_methylation_expression(
    tumor = "LUAD", 
    gene = "xxx"), 
    "No valid")
  
  ## returns an error if more than one gene entered  
  expect_error(prepare_TCGA_methylation_expression(
    tumor = "LUAD", 
    gene = c("MAGEA1", "MAGEA3")), 
    "No valid")
  
  ## Collects the right methylation values 
  res <- prepare_TCGA_methylation_expression(tumor = "ESCA", 
                                             gene = "MAGEA1")
  # expected_res <- dput(res)
  # save(expected_res, file = "tests/testthat/methylation_expression_MAGEA1_in_ESCA.rda")
  load("tests/testthat/methylation_expression_MAGEA1_in_ESCA.rda")
  expect_identical(res, expected_res)
  
  ## Only tumor samples are returned when include_normal_tissues is FALSE
  res <- prepare_TCGA_methylation_expression(tumor = "ESCA", 
                                             gene = "MAGEA1",
                                             include_normal_tissues = FALSE)
  
  expect_true(!"Peritumoral" %in% res$tissue)
  
  ## Peritumoral samples are included when include_normal_tissues is TRUE
  res <- prepare_TCGA_methylation_expression(tumor = "ESCA", 
                                             gene = "MAGEA1",
                                             include_normal_tissues = TRUE)
  
  expect_true("Peritumoral" %in% res$tissue)
  
  
  
})
