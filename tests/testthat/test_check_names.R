test_that("check_names() works", {
    input1 <- "var1"
    valid1 <- c("var1", "var2")
    res1 <- "var1"
    expect_identical(CTexploreR:::check_names(input1, valid1),
                     res1)
    input2 <- c("var0", "var1")
    res2 <- "var1"
    expect_identical(CTexploreR:::check_names(input2, valid1),
                     res2)
    expect_warning(CTexploreR:::check_names(input2, valid1),
                   "names invalid")
    input3 <- c("var2", "var1")
    valid3 <- paste0("var", 1:10)
    ## order is preserved
    expect_equal(CTexploreR:::check_names(input3, valid3),
                 c("var2", "var1"))
    expect_false(identical(CTexploreR:::check_names(input3, valid3),
                           c("var1", "var2")))
    ## empty result
    expect_identical(CTexploreR:::check_names("var", valid3),
                     character())
    expect_warning(CTexploreR:::check_names("var", valid3),
                     "names invalid")
})
