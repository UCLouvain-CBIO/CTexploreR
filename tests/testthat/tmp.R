test_that("GTEX_expression() returns the correct number gene", {
    ## 2 genes gives a matrix of 2 rows
    my_genes <- c("MAGEA1", "MAGEA3")
    return_val <- TRUE
    res <- CTexploreR:::GTEX_expression(my_genes, return = return_val)
    
    expect_identical(nrow(res), 2L)
    expect_equal(nrow(res), 2) # idem mais pas besoin du 2L
    #expect_identical(rownames(res), my_genes)
    expect_identical(sort(rownames(res)), sort(my_genes))
    #expect_true(res, "matrix")
    expect_true(inherits(res, "matrix"))
    expect_type(res, "double")
    # dput(res) à copier .... hehehe
    res0 <- structure(
      c(12.8107, 10.6212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0), 
      dim = c(2L, 32L), 
      dimnames = list(
        c("MAGEA3", "MAGEA1"), c("Testis", "Ovary", "Adipose", "Adrenal Gland",
                                 "Artery", "Bladder", "Blood", "Brain", "Breast", "Cervix",
                                 "Colon", "Esophagus", "Fallopian Tube", "Fibroblasts", "Heart",
                                 "Kidney", "Liver", "Lung", "Lymphocytes_EBV", "Muscle", "Nerve",
                                 "Pancreas", "Pituitary", "Prostate", "Salivary Gland", "Skin",
                                 "Small Intestine", "Spleen", "Stomach", "Thyroid", "Uterus",
                                 "Vagina")))
    expect_identical(res, res0)
})

