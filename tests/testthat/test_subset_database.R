test_that("subset_database() works", {
  dataset <- CTdata::GTEX_data()
  var <- c("MAGEA1", "MAGEA3")
  res <- CTexploreR:::subset_database(var, dataset)
  
  # Check two genes gives two genes in dataset
  expect_identical(sort(rowData(res)$external_gene_name),
                   sort(var))
  expect_equal(nrow(res), 2)
  
  # Check no variable gives CT_genes
  res_no <- CTexploreR:::subset_database(data = dataset)
  expect_identical(rowData(res_no)$external_gene_name,
                   CTdata::CT_genes()$external_gene_name)
  expect_equal(nrow(res_no), nrow(CTdata::CT_genes()))

  
  # Check when fake gene name gives nothing + warning
  var_fake <- c("not_existing")
  expect_warning(CTexploreR:::subset_database(var_fake, dataset))
  res_fake <- CTexploreR:::subset_database(var_fake, dataset)
  expect_equal(nrow(res_fake), 0)
})
