test_that("normal_tissue_expression_multimapping() works", {
  
  ## returns a matrix of double when Return is set to TRUE
  ## returns a warning when an invalid gene is entered
  ## n valid genes in input returns a matrix of n expected rownames
  ## Test multimapping = TRUE
  tested_genes <- c("MAGEA3", "MAGEA6", "")
  expect_warning(res <- normal_tissue_expression_multimapping(
    genes = tested_genes, 
    multimapping = TRUE, return = TRUE), "names invalid")
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 2) 
  tested_genes <- c("MAGEA3", "MAGEA6")
  ensembl_ids <- c("ENSG00000221867", "ENSG00000197172")
  expected_assay <- assay(normal_tissues_multimapping_data(), 
                          "TPM_with_multimapping")[ensembl_ids, ]
  rownames(expected_assay) <- tested_genes
  expect_equal(res[tested_genes,], expected_assay[tested_genes,])
  
  ## Test the "log_TPM" units argument
  res_in_log <- normal_tissue_expression_multimapping(
    genes = tested_genes, multimapping = TRUE, 
    units = "log_TPM", return = TRUE)
  expect_equal(res_in_log, log1p(res)) 
  
  ## Test multimapping = FALSE
  res <- normal_tissue_expression_multimapping(
    genes = tested_genes, multimapping = FALSE, return = TRUE)
  expected_assay <- assay(normal_tissues_multimapping_data(), 
                          "TPM_no_multimapping")[ensembl_ids, ]
  rownames(expected_assay) <- tested_genes
  expect_equal(res[tested_genes,], expected_assay[tested_genes,])
  
  ## Test that the function returns a heatmap by default
  ## res <- normal_tissue_expression_multimapping(
  ##  genes = c("MAGEA3", "MAGEA6"), multimapping = TRUE)
  ## expect_s4_class(res, "Heatmap")
  ## vdiffr::expect_doppelganger("normal_tissue_expression_multimapping_on_MAGE", 
  ##                            fig = res)
})

