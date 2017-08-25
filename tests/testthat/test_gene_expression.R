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

test_that("getTxDensity returns list", {
    expect_type(getTxDensity(tx, annox), "list")
})

test_that("getTxDensity joins broken transcript", {
    result <- getTxDensity(tx, annox)
    names(result) <- NULL
    expected <- list(
        0,                              # No head signal
        1,                              # Perfect overlap after join
        1,                              # Tail signal (why is it not 0?)
        1                               # Perfect overlap after join
    )
    expect_equal(result, expected)
})

test_that("getTxDensity splits combined transcript", {
    result <- head(getTxDensity(annox, tx), 3)
    names(result) <- NULL
    expected <- list(
        0,                                 # No head signal
        head(width(annox) / width(tx), 1), # Overlap % of first transcript
        0                                  # No tail signal
    )
    expect_equal(result, expected, tolerance = 1e-3)
})

test_that("getLIValues linearly interpolates single value", {
    vals <- 9
    n <- 5
    expect_equal(getLIValues(vals, n), rep(vals, n))
})

test_that("evaluateHMMInAnnotations returns data.frame", {
    expect_type(evaluateHMMInAnnotations(tx, annox), "list")
})

test_that("evaluateHMMInAnnotations reports merged", {
    result <- evaluateHMMInAnnotations(tx, annox)$eval
    expect_true(result$merged == 1)
})

test_that("evaluateHMMInAnnotations reports dissociated", {
    result <- evaluateHMMInAnnotations(annox, tx)$eval
    expect_true(result$dissociated == 1)
})
