test_that("subset_database() works", {
  dataset <- CTdata::GTEX_data()
  var <- c("MAGEA1", "MAGEA3")
  
  ## Returns a Summarized Experiment
  res <- CTexploreR:::subset_database(var, dataset)
  expect_s4_class(res, "SummarizedExperiment")
  
  ## Check two genes gives two genes in dataset
  res <- CTexploreR:::subset_database(var, dataset)
  expect_identical(sort(rowData(res)$external_gene_name),
                   sort(var))
  expect_equal(nrow(res), 2)
  
  ## Check no variable gives strict CT_genes
  res_no <- CTexploreR:::subset_database(data = dataset)
  expect_identical(rowData(res_no)$external_gene_name,
                   CT_genes[CT_genes$CT_gene_type == "CT_gene", 
                            "external_gene_name", drop = TRUE])
  expect_equal(nrow(res_no), nrow(CT_genes[CT_genes$CT_gene_type == "CT_gene", ]))
  
  ## Check when fake gene name gives nothing + warning
  var_fake <- c("not_existing")
  expect_warning(res_fake <- CTexploreR:::subset_database(var_fake, dataset))
  expect_equal(nrow(res_fake), 0)
})
