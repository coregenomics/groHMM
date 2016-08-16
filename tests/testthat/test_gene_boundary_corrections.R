library(groHMM)
context("Gene boundary corrections")

## Fixtures
tx <- GRanges(c("chr7:1000-30000"), strand="+")
annox <- GRanges(c("chr7:1000-11000",
                   "chr7:20000-30000"), strand="+")

## Introspect default values
gap <- formals(breakTranscriptsOnGenes)$gap

test_that("breakTranscriptsOnGenes returns (gap - 1) sized IRange", {
    expect_equal(
        ranges(
            breakTranscriptsOnGenes(tx, annox)[1]),
        narrow(
            ranges(
                GRanges(c("chr7:1000-20000"))),
            end=-gap - 1))
})
