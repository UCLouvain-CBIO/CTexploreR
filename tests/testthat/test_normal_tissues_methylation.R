test_that("normal_tissues_methylation() works", {
    
    ## returns a matrix of double
    res <- normal_tissues_methylation(gene = "TDRD1", 1000, 0, return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## the return matrix has the expected number of rows
    expect_equal(nrow(res), ncol(CT_methylation_in_tissues())) 
    
    ## No valid gene name returns an error
    expect_error(normal_tissues_methylation(gene = "xxx", 1000, 0), 
                 "is not in the CT database")
    
    ## a warning message specifies that the max nt_up / nt_down must be < 5000
    expect_warning(normal_tissues_methylation(gene = "TDRD1", 
                                              nt_up = 20000, nt_down = 20000),
                   "maximum number")
    
    ## Collects the right methylation values (for a gene transcribed in sense)
    res <- normal_tissues_methylation(gene = "TDRD1", 1000, 200, return = TRUE)
    TSS <- 114179353
    promoter_gr <- GRanges(seqnames = "chr10",
                           strand = "+",
                           ranges = IRanges(
                             start = TSS - 1000,
                             end = TSS + 200))
    exp_data <- subsetByOverlaps(CT_methylation_in_tissues(), promoter_gr)
    exp_res <- as.matrix(assay(exp_data))
    rownames(exp_res) <- exp_data@rowRanges@ranges@start -TSS
    exp_res <- t(exp_res)
    expect_identical(res, exp_res)
    
    ## Collects the right methylation values (for a gene transcribed in antisense)
    res <- normal_tissues_methylation(gene = "SSX3", 1000, 200, return = TRUE)
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

    ## The function works even when no CpG whithin the selected range
    res <- normal_tissues_methylation(gene = "MAGEA1", 5, 5, return = TRUE)
    expect_equal(nrow(res), 0) 
    
    ## Test that the function returns a heatmap
    res <- normal_tissues_methylation(gene = "MAGEA1")
    expect_s4_class(res, "Heatmap")
})
