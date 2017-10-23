context("Gene boundary corrections")

## Fixtures
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
                        end = -gap - 1))
})

test_that("breakTranscriptsOnGenes honors negative strand", {
    strand_ <- "-"
    tx_ <- tx
    annox_ <- annox
    strand(tx_) <- strand_
    strand(annox_) <- strand_
    expect_equal(ranges(breakTranscriptsOnGenes(tx_, annox_,
                                                strand = strand_)[1]),
                 narrow(ranges(GRanges(c("chr7:1000-11000")))))
})

test_that("breakTranscriptsOnGenes honors no annotation breaks", {
    annox_ <- shift(tx, width(tx))
    expect_equal(ranges(breakTranscriptsOnGenes(tx, annox_)), ranges(tx))
})

test_that("breakTranscriptsOnGenes honors multiple annotation breaks", {
    annox_ <- GRanges(c("chr7:1000-11000",
                        "chr7:12000-20000",
                        "chr7:20001-30000"), strand="+")
    expect_equal(length(breakTranscriptsOnGenes(tx, annox_)), 3)
})

test_that("breakTranscriptsOnGenes honors multiple transcripts", {
    tx_ <- GRanges(c("chr7:1000-12000",
                     "chr7:20000-30000"),
                   strand="+")
    annox_ <- GRanges(c("chr7:1001-6000",
                        "chr7:6005-12000",
                        "chr7:20001-25000",
                        "chr7:25005-30000"),
                      strand="+")
    expect_equal(length(breakTranscriptsOnGenes(tx_, annox_,
                                                geneSize = 4995)), 4)
})

test_that("Overlapping annotations raise error", {
    tx_ <- GRanges("chr1:38834900-38873499:-")
    annox_ <- GRanges(c("chr1:38838198-38874494:-",
                        "chr1:38864236-38873368:-"))
    expect_error(breakTranscriptsOnGenes(tx_, annox_, strand = "-"),
                 "overlap")
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

    expect_equal(numTxCb(txShortest), 1)
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

test_that("breakInterval always returns an S4 GRanges object when bounded", {
    ## Empty objects
    expect_s4_class(breakInterval(GRanges(), NULL), "GRanges")
    expect_s4_class(breakInterval(tx, NULL), "GRanges")
    expect_s4_class(breakInterval(tx, start(annox[2])), "GRanges")
    expect_equal(length(breakInterval(GRanges(), NULL)), 0)
    expect_equal(length(breakInterval(tx, NULL)), 1)
    expect_equal(length(breakInterval(tx, start(annox[2]))), 2)
})

test_that("breakInterval errors when 'brPos' is out of 'gr' Range", {
    expect_error(breakInterval(tx, start(tx) - 500))
    expect_error(breakInterval(tx, end(tx) + 500))
})

test_that("breakInterval honors 'strand'", {
    gap <- formals(breakInterval)$gap
    brPos <- start(annox[2])

    res <- breakInterval(tx, brPos, strand="+")
    expect_equal(end(res)[1], brPos - gap)
    expect_equal(start(res)[2], brPos)

    res <- breakInterval(tx, brPos, strand="-")
    expect_equal(end(res)[1], brPos)
    expect_equal(start(res)[2], brPos + gap)
})

test_that("makeConsensusAnnotations errors on missing 'gene ID'", {
    expect_error(makeConsensusAnnotations(tx), "Missing gene .* column")
})

test_that("makeConsensusAnnotations output is non-overlapping and condensed", {
    annox_ <- GRanges(
        c(
            "chr7:1000-6000:+",         # multi overlapping isoform
            "chr7:2000-6000:+",         # multi overlapping isoform
            "chr7:1000-12000:+",        # multi overlapping isoform
            "chr7:15000-17000:+",       # single overlapping isoform
            "chr7:16000-17000:+",       # single overlapping isoform
            "chr7:20000-22000:+",       # non-overlapping isoform
            "chr7:22010-23000:+",       # non-overlapping isoform
            "chr7:25000-26000:+",       # mixed strand isoform
            "chr7:25000-27000:-",       # mixed strand isoform
            "chr7:29000-32000:+",       # mixed chrom isoform
            "chr8:29000-33000:+",       # mixed chrom isoform
            "chr7:40000-45000:+",       # + overlapping annotation
            "chr7:44000-50000:+",       # + overlapping annotation
            "chr7:50000-55000:-",       # - overlapping annotation
            "chr7:54000-60000:-"        # - overlapping annotation
        ),
        gene_id=c(rep(1, 3), rep(2:5, each=2), 6:9))
    result <- sort(makeConsensusAnnotations(annox_))
    names(result) <- NULL
    expected <- sort(GRanges(
        c(
            "chr7:1000-6000:+",         # most hit isoform for multi overlapping
            "chr7:16000-17000:+",       # shorter isoform for single overlapping
            "chr7:20000-22000:+",       # longer isoform for non-overlapping
            "chr7:25000-27000:-",       # longer isoform for mixed strand
            "chr8:29000-33000:+",       # longer isoform for mixed chrom
            "chr7:40000-43999:+",       # trimmed + overlapping annotation
            "chr7:44000-50000:+",       # unchanged + overlapping annotation
            "chr7:50000-55000:-",       # unchanged - overlapping annotation
            "chr7:55001-60000:-"        # trimmed - overlapping annotation
        ),
        gene_id=1:9))
    expect_s4_class(result, "GRanges")
    expect_true(isDisjoint(result))
    expect_equal(result, expected)
})

test_that("makeConsensusAnnotations drops missing gene ID values", {
    annox_ <- GRanges(c("chr7:1000-2000", "chr7:3000-4000"),
                      gene_id=CharacterList(list(character(0), "1")))
    expect_warning(makeConsensusAnnotations(annox_), "dropped")
})

test_that("makeConsensusAnnotations drops multiple gene ID values", {
    annox_ <- GRanges(c("chr7:1000-2000", "chr7:3000-4000"),
                      gene_id=CharacterList(list(c("1", "2"), "1")))
    expect_warning(makeConsensusAnnotations(annox_), "dropped")
})
