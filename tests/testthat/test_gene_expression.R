context("Gene expression")

test_that("expressedGenes returns data.frame", {
    expect_s3_class(expressedGenes(annox, tx), "data.frame")
    expect_warning(expressedGenes(annox, GRanges()), "empty")
})

test_that("expressedGenes validates reads", {
    expect_error(expressedGenes(annox, "reads"), "GRanges")
    expect_error(expressedGenes(annox, 1:10), "GRanges")
    expect_error(expressedGenes(annox, GRangesList()), "GRanges")
})

test_that("expressedGenes propagates mcols(features) to output data.frame", {
    expected <- colnames(mcols(annox))
    result <- colnames(expressedGenes(annox, tx))
    cols_common <- intersect(result, expected)
    expect_equal(expected, cols_common)
})
