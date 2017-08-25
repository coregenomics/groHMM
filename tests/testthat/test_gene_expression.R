context("Gene expression")

reads_ <- GRanges("chr7",
                  IRanges(start=c(1000:1003, 1100:1101),
                          width=rep(1, 6)),
                  strand=rep(c("+", "-"), 3))
features <- GRanges(c("chr7:1000-1000:+",
                      "chr7:1001-1001:-"))
values <- c(0, 0.5, 1.5, 0.5, 0)
expected <- list(sense = Rle(values, c(7, 2, 4, 2, 5)),
                 antisense = Rle(values, c(6, 2, 4, 2, 6)))


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

test_that("pausingIndex handles empty inputs", {
    gr <- tx

    result <- pausingIndex(GRanges(), GRanges())
    expect_s3_class(result, "data.frame")
    expect_true(NROW(result) == 0)

    result <- pausingIndex(GRanges(), gr)
    expect_s3_class(result, "data.frame")
    expect_true(NROW(result) == 0)

    result <- pausingIndex(gr, GRanges())
    expect_s3_class(result, "data.frame")
    expect_true(NROW(result) == 1)
    result$GeneID <- NULL
    expect_true(all(result == 0))
})

test_that("pausingIndex returns correct result ", {
    ## Use man page example
    features <- GRanges("chr7:2394474-2420377:+")
    result <- pausingIndex(features, reads[[1]])
    expected <- data.frame(
        Pause = 0.5,
        Body = 0.01260892,
        Fisher = 5.257256e-7,
        GeneID = factor("1"),
        CIlower = 25.87514,
        CIupper = 60.68879,
        PauseCounts = 25,
        BodyCounts = 313,
        uPCounts = 1,
        uGCounts = 337
    )
    expect_equal(result, expected, tolerance = 1e-6)

    ## Labelled genes
    symbol <- "my_gene"
    mcols(features) <- DataFrame(symbol = symbol)
    result <- pausingIndex(features, reads[[1]])
    expected$GeneID <- factor(symbol)
    expect_equal(result, expected, tolerance = 1e-6)
})

test_that("runMetaGene validates anchorType", {
    expect_error(runMetaGene(annox, tx, anchorType="bad_value"), "anchorType")
})

test_that("runMetaGene returns correct run length encodings", {
    ## Use example
    result <- runMetaGene(features, reads_, size=4, up=10)
    expect_equal(result, expected)

    result <- runMetaGene(features, reads_, size=4, up=10, anchorType="TTS")
    expect_equal(result, expected)
})

test_that("metaGene ensures reads or coverage as input", {
    expect_error(metaGene(annox, reads=NULL), "reads")
    expect_error(metaGene(annox, plusCVG=NULL), "plusCVG")
    expect_error(metaGene(annox, minusCVG=NULL), "minusCVG")
})

test_that("metaGene accepts 'reads' combined input", {
    expect_s4_class(metaGene(features, reads_, size=4, up=10), "Rle")
})
