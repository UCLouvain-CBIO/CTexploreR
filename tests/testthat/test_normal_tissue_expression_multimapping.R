test_that("normal_tissue_expression_multimapping() works", {
    
    ## returns a matrix of double when Return is set to TRUE
    res <- normal_tissue_expression_multimapping(
      genes = c("MAGEA3"), 
      multimapping = TRUE, return = TRUE)
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    
    ## Test that the function returns a heatmap by default
    res <- normal_tissue_expression_multimapping(
      genes = c("MAGEA3", "MAGEA6"), multimapping = TRUE)
    expect_s4_class(res, "Heatmap")
    
    ## n valid genes in input returns a matrix of n expected rownames
    res <- normal_tissue_expression_multimapping(
      genes = c("MAGEA3", "MAGEA6"), 
      multimapping = TRUE, return = TRUE)
    expect_equal(nrow(res), 2) 
    expect_identical(sort(rownames(res)), sort(c("MAGEA3", "MAGEA6")))
    
    ## works with only one gene in input
    res <- normal_tissue_expression_multimapping(
      genes = c("MAGEA3"), 
      multimapping = TRUE, return = TRUE)
    expect_equal(nrow(res), 1) 
    
    ## works but returns a warning when no valid gene is entered
    res <- normal_tissue_expression_multimapping(
      genes = "", multimapping = TRUE, return = TRUE)
    expect_equal(nrow(res), 0)
    expect_warning(normal_tissue_expression_multimapping(
      genes = "", multimapping = TRUE, return = TRUE), "names invalid")
    
    ## collects correctly data not multimapped
    tested_genes <- c("MAGEA3", "MAGEA6")
    res <- normal_tissue_expression_multimapping(
      genes = tested_genes, multimapping = FALSE, return = TRUE)
    CT_ensembl_gene <- as.data.frame(
      CT_genes[CT_genes$external_gene_name %in% tested_genes, ])
    rownames(CT_ensembl_gene) <- CT_ensembl_gene$external_gene_name
    asked_ensembl_ids <- CT_ensembl_gene[tested_genes, ]$ensembl_gene_id
    expected_assay <- assay(normal_tissues_multimapping_data(), 
                            "TPM_no_multimapping")[asked_ensembl_ids, ]
    rownames(expected_assay) <- tested_genes
    expect_equal(res[tested_genes,], expected_assay[tested_genes,])
    
    ## Tests if collects correctly multimapping data
    res <- normal_tissue_expression_multimapping(
      genes = tested_genes, multimapping = TRUE, return = TRUE)
    expected_assay <- assay(normal_tissues_multimapping_data(), 
                            "TPM_with_multimapping")[asked_ensembl_ids, ]
    rownames(expected_assay) <- tested_genes
    expect_equal(res[tested_genes,], expected_assay[tested_genes,])
    
    ## No valid multimapping paramater returns an error
    expect_error(normal_tissue_expression_multimapping(
      genes = tested_genes, return = TRUE), "multimapping parameter")

    ## Test the "log_TPM" units argument
    res_in_TPM <- normal_tissue_expression_multimapping(
      genes = tested_genes, multimapping = TRUE, return = TRUE)
    res_in_log <- normal_tissue_expression_multimapping(
      genes = tested_genes, multimapping = TRUE, units = "log_TPM", 
      return = TRUE)  
    expect_equal(res_in_log, log1p(res_in_TPM)) 
    
})
