test_that("HPA_cell_type_expression() works", {

  ## returns a SingleCellExperiment when Return is set to TRUE
  ## n valid genes in input returns a matrix of n expected rownames
  ## rownames in the returned data correspond to ENSEMBL ids
  expect_warning(res <- HPA_cell_type_expression(c("MAGEA1", "MAGEA3", "xxx"),
                                                 return = TRUE),
                 "valid types")
  expect_s4_class(res, "SingleCellExperiment")
  expect_equal(nrow(res), 2)
  expect_true(all(unique(rownames(res)) %in%
                    c("ENSG00000198681", "ENSG00000221867")))

  ## Test that the function returns a heatmap by default
  ## works with only one gene in input
  ## res <- HPA_cell_type_expression("MAGEA1")
  ## expect_s4_class(res, "Heatmap")
  ## vdiffr::expect_doppelganger("HPA_cell_type_expression_MAGEA1", fig = res)
})
