## Boilerplate for enforcing lintr with testthat from lintr README.md
if (requireNamespace("lintr", quietly = TRUE)) {
    context("Code style")
    test_that("Package Style", {
        lintr::expect_lint_free()
    })
}
