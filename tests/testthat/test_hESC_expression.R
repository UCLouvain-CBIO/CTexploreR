test_that("hESC_expression() works", {
  
  ## returns a SummarizedExperiment
  res <- hESC_expression(genes = c("MAGEA1", "MAGEA3"), 
                           values_only = TRUE)
  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(nrow(res), 2) 
  expect_identical(sort(rowData(res)$external_gene_name), 
                   sort(c("MAGEA1", "MAGEA3")))
  
  ## Test the "log_TPM" units argument
  res_in_log <- hESC_expression(c("MAGEA1", "MAGEA3"), 
                                units = "log_TPM", 
                                values_only = TRUE)
  expect_equal(assay(res_in_log)[, "H1"], log1p(assay(res)[, "H1"])) 
  
  # Check default genes with CTP
  res_includeCTP <- hESC_expression(include_CTP = TRUE, values_only = TRUE)
  expect_equal(nrow(res_includeCTP), nrow(CT_genes))
  
  ## Test that the function returns a heatmap
  ## res <- hESC_expression(genes = c("MAGEA1", "MAGEA3"))
  ## expect_s4_class(res, "Heatmap")
})