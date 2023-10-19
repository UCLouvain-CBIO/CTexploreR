test_that("normal_tissues_methylation() works", {
  
  
  ## returns a matrix of double when values_only is set to TRUE
  res <- normal_tissues_methylation(gene = "TDRD1", 500, 0, values_only = TRUE)
  expect_true(inherits(res, "matrix"))
  expect_type(res, "double")
  expect_equal(nrow(res), 14) 
  
  ## Collects the right methylation values (for a gene transcribed in sense)
  ## tests the nt_up and nt_down parameters
  TSS <- 114179353
  promoter_gr <- GRanges(seqnames = "chr10",
                         strand = "+",
                         ranges = IRanges(
                           start = TSS - 500,
                           end = TSS - 0))
  exp_data <- subsetByOverlaps(CT_methylation_in_tissues(), promoter_gr)
  exp_res <- as.matrix(assay(exp_data))
  rownames(exp_res) <- exp_data@rowRanges@ranges@start -TSS
  exp_res <- t(exp_res)
  expect_identical(res, exp_res)
  
  ## Collects the right methylation values (for a gene transcribed in antisense)
  res <- normal_tissues_methylation(gene = "SSX3", values_only = TRUE)
  TSS <- 48356703
  promoter_gr <- GRanges(seqnames = "chrX",
                         strand = "-",
                         ranges = IRanges(
                           start = TSS - 200,
                           end = TSS + 1000))
  exp_data <- subsetByOverlaps(CT_methylation_in_tissues(), promoter_gr)
  exp_res <- as.matrix(assay(exp_data))
  rownames(exp_res) <- TSS - exp_data@rowRanges@ranges@start
  exp_res <- t(exp_res)
  expect_identical(res, exp_res)
  
  ## Works even when no CpG whithin the selected range
  res <- normal_tissues_methylation(gene = "MAGEA1", 5, 5, values_only = TRUE)
  expect_equal(nrow(res), 0) 
  
  ## returns a heatmap by default
  ## a warning message specifies that the max nt_up / nt_down must be < 5000
  ## expect_warning(res <- normal_tissues_methylation(gene = "TDRD1", 
  ##                                           nt_up = 20000, nt_down = 20000),
  ##                "maximum number")
  ## expect_s4_class(res, "Heatmap")
  ## vdiffr::expect_doppelganger("normal_tissues_methylation_on_TDRD1", fig = res)
})
