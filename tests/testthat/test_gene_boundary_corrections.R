library(groHMM)
context("Gene boundary corrections")

## Fixtures
tx <- GRanges(c("chr7:1000-30000"), strand="+")
annox <- GRanges(c("chr7:1000-11000",
                   "chr7:20000-30000"), strand="+")
txC <- annox
annoxC <- tx
numTxBr <- function(annox) length(ranges(breakTranscriptsOnGenes(tx, annox)))
numTxCb <- function(txC) length(ranges(combineTranscripts(txC, annoxC)))

test_that("Inputs 'tx' and 'annox' are validated", {
    expect_error(breakTranscriptsOnGenes(NULL, annox),
                 "'tx'", ignore.case = TRUE)
    expect_error(breakTranscriptsOnGenes(GRanges(), annox),
                 "'tx'", ignore.case = TRUE)
    expect_error(combineTranscripts(NULL, annox),
                 "'tx'", ignore.case = TRUE)
    expect_error(combineTranscripts(GRanges(), annox),
                 "'tx'", ignore.case = TRUE)

    expect_error(breakTranscriptsOnGenes(tx, NULL),
                 "'annox'", ignore.case = TRUE)
    expect_warning(breakTranscriptsOnGenes(tx, GRanges()),
                   "'annox'", ignore.case = TRUE)
    expect_error(combineTranscripts(tx, NULL),
                 "'annox'", ignore.case = TRUE)
    expect_warning(combineTranscripts(tx, GRanges()),
                   "'annox'", ignore.case = TRUE)
})

test_that("breakTranscriptsOnGenes returns 'gap' sized IRanges", {
    gap <- formals(breakTranscriptsOnGenes)$gap

    expect_equal(ranges(breakTranscriptsOnGenes(tx, annox)[1]),
                 narrow(ranges(GRanges(c("chr7:1000-20000"))),
                        end=-gap - 1))
})

test_that("Gene repair functions honor 'geneSize'", {
    geneSize <- formals(breakTranscriptsOnGenes)$geneSize
    annoxShortest <- resize(annox, geneSize + 1)
    annoxTooShort <- narrow(annoxShortest, end=-2)

    expect_equal(numTxBr(annoxShortest), 2)
    expect_equal(numTxBr(annoxTooShort), 1)

    geneSize <- formals(combineTranscripts)$geneSize
    txShortest <- resize(txC, geneSize + 1)
    txTooShort <- narrow(txShortest, end=-2)

    expect_equal(numTxCb(txShortest), 2)
    expect_equal(numTxCb(txTooShort), 1)
})

test_that("breakTranscriptsOnGenes honors 'threshold'", {
    threshold <- formals(breakTranscriptsOnGenes)$threshold
    annoxMinOverlap <- shift(annox, width(annox) * (1 - threshold))
    annoxInsuffOverlap <- shift(annox, width(annox) * (1.05 - threshold))

    expect_equal(numTxBr(annoxMinOverlap), 2)
    ## The second gene now does not overlap with the transcript.
    expect_equal(numTxBr(annoxInsuffOverlap), 1)
})
