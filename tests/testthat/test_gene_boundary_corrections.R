library(groHMM)
context("Gene boundary corrections")

## Fixtures
tx <- GRanges(c("chr7:1000-30000"), strand="+")
annox <- GRanges(c("chr7:1000-11000",
                   "chr7:20000-30000"), strand="+")
numTx <- function(annox) length(ranges(breakTranscriptsOnGenes(tx, annox)))

test_that("breakTranscriptsOnGenes checks validity of 'tx'", {
    expect_error(breakTranscriptsOnGenes(NULL, annox),
                 "'tx'", ignore.case = TRUE)
    expect_error(breakTranscriptsOnGenes(GRanges(), annox),
                 "'tx'", ignore.case = TRUE)
})

test_that("breakTranscriptsOnGenes checks validity of 'annox'", {
    expect_error(breakTranscriptsOnGenes(tx, NULL),
                 "'annox'", ignore.case = TRUE)
    expect_warning(breakTranscriptsOnGenes(tx, GRanges()),
                   "'annox'", ignore.case = TRUE)
})

test_that("breakTranscriptsOnGenes returns 'gap' sized IRanges", {
    gap <- formals(breakTranscriptsOnGenes)$gap

    expect_equal(ranges(breakTranscriptsOnGenes(tx, annox)[1]),
                 narrow(ranges(GRanges(c("chr7:1000-20000"))),
                        end=-gap - 1))
})

test_that("breakTranscriptsOnGenes honors 'geneSize'", {
    geneSize <- formals(breakTranscriptsOnGenes)$geneSize
    annoxShortest <- resize(annox, geneSize + 1)
    annoxTooShort <- narrow(annoxShortest, end=-2)

    ## Minimum acceptable annotation geneSize should work just as well
    ## as a longer annotation.
    expect_equal(numTx(annoxShortest), 2)
    ## Annotation should be ignored when geneSize is too short.
    expect_equal(numTx(annoxTooShort), 1)
})

test_that("breakTranscriptsOnGenes honors 'threshold'", {
    threshold <- formals(breakTranscriptsOnGenes)$threshold
    annoxMinOverlap <- shift(annox, width(annox) * (1 - threshold))
    annoxInsuffOverlap <- shift(annox, width(annox) * (1.05 - threshold))

    expect_equal(numTx(annoxMinOverlap), 2)
    ## The second gene now does not overlap with the transcript.
    expect_equal(numTx(annoxInsuffOverlap), 1)
})
