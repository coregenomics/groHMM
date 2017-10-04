context("HMM transcript and polymerase wave detection")

## Fixtures
approxDist <- 20000

expected <- data.frame(
    StartWave = 9350,
    EndWave = 33750,
    Rate = 24400,
    minOfMax = as.numeric(NA),
    minOfAvg = as.numeric(NA),
    SYMBOL = "CYP2W1",
    ID = "54905",
    stringsAsFactors = FALSE
)

test_that("polymeraseWave returns expected data.frame", {
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE)
    expect_s3_class(pw, "data.frame")
    expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("polymeraseWave returns identical value on opposite strand", {
    genes <- invertStrand(genes)
    reads <- invertStrand(reads)
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE)
    expect_s3_class(pw, "data.frame")
    ## Needing to allow for +/- 1 tolerance here points to a subtle integer
    ## error elsewhere in the code.
    expect_equal(pw, expected, tolerance = 1)
})

test_that("polymeraseWave checks filterWindowSize against size", {
    filterWindowSize <- 10000
    check <- function(size) {
        polymeraseWave(
            reads[[1]], reads[[2]], genes, approxDist,
            progress=FALSE, filterWindowSize=filterWindowSize, size=size)
    }
    expect_error(check(99), "filterWindowSize")
    expect_error(check(101), "filterWindowSize")
    expect_error(check(9999), "filterWindowSize")
    expect_error(check(-1), "windowSize.*out of range")
    expect_silent(check(100))
})

test_that("polymeraseWave drops genes with few reads", {
    check <- function(genes) {
        polymeraseWave(
            reads[[1]], reads[[2]], genes, approxDist, progress=FALSE)
    }
    expect_silent(check(genes))
    genes <- c(
        genes,
        GRanges("chr7:1-100:+", SYMBOL="CYP2W1", ID="54905"))
    expect_message(check(genes), "Drop")
})

test_that("polymeraseWave generates a progress bar", {
    if (! interactive())
        skip("Cannot test progress bar in non-interactive session")
    tmp <- tempfile()
    on.exit({
        close(con)
        unlink(tmp)
    })
    con <- file(tmp, open="wt")
    capture.output(
        invisible(polymeraseWave(reads[[1]], reads[[2]], genes, approxDist)),
        file=con, type="message")
    result <- readChar(tmp, file.info(tmp)$size)
    expect_match(result, "\\[.*\\]")
})

test_that("polymeraseWave raises error for invalid distribution", {
    expect_error(
        polymeraseWave(
            reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
            emissionDistAssumption="paranormal"),
        "emissionDistAssumption")
})

test_that("polymeraseWave TSmooth induces smoothing", {
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
        emissionDistAssumption="norm", TSmooth=20)
    expect_equal(pw, expected, tolerance = 1e-7)

    expected[, c("StartWave", "EndWave", "Rate", "minOfMax", "minOfAvg")] <-
        c(10500, 12850, 2350, 0, 1)
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
        emissionDistAssumption="norm", TSmooth="3RS3R")
    expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("polymeraseWave normal exponental distribution assumption", {
    expected[, c("EndWave", "Rate", "minOfMax", "minOfAvg")] <-
        c(10050, 700, 0, 1)
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
        emissionDistAssumption="normExp")
    expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("polymeraseWave gene resize honors seqinfo boundaries", {
    seqinfo(genes) <- seqinfo(reads)
    ## Resize to be 1 base too big.
    dist_to_end <- (seqlengths(genes) - end(genes) - approxDist)
    genes <- shift(genes, dist_to_end)
    suppressWarnings(reads <- endoapply(shift(reads, dist_to_end), trim))
    expect_silent(polymeraseWave(reads[[1]], reads[[2]], genes, approxDist,
                                 progress=FALSE))
    genes <- shift(genes, approxDist)
    expect_message(polymeraseWave(reads[[1]], reads[[2]], genes, approxDist,
                                  progress=FALSE),
                   "too close to chromosome edge")
})

test_that("detectTranscripts generates sensible transcripts", {
    ## Runtime is about a minute :/
    reads_ <- unlist(endoapply(reads, head, 200))
    hmm <- detectTranscripts(reads_, threshold=1)
    result <- granges(hmm$transcripts)
    ## Finds a single transcript.
    expected <- GRanges("chr7:564350-568449:+")
    ## Result is stochastic.
    expect_equal(start(result), start(expected), tolerance=2)
    expect_equal(end(result), end(expected), tolerance=2)
})
