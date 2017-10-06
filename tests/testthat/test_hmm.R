context("HMM transcript and polymerase wave detection")

## Fixtures
approxDist <- 20000
expected <- list(GRanges(
    seqnames = "chr7",
    ranges = IRanges(
        2394474 - 9350,
        2420377 - 10853),
    strand = "+",
    SYMBOL = "CYP2W1",
    ID = "54905",
    rate = as.numeric(24400),
    min_of_max = as.numeric(NA),
    min_of_avg = as.numeric(NA)
))

test_that("polymeraseWave returns expected GRanges", {
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE)
    expect_s4_class(pw[[1]], "GRanges")
    expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("polymeraseWave returns identical value on opposite strand", {
    genes <- invertStrand(genes)
    reads <- invertStrand(reads)
    expected[[1]] <- invertStrand(expected[[1]])
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE)
    expect_s4_class(pw[[1]], "GRanges")
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
    bpparam <- SerialParam(stop.on.error = FALSE)
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
        emissionDistAssumption="paranormal", BPPARAM=bpparam)
    expect_false(bpok(pw))
    expect_match(pw[[1]]$message, "emissionDistAssumption")
})

test_that("polymeraseWave TSmooth induces smoothing", {
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
        emissionDistAssumption="norm", TSmooth=20)
    expect_equal(pw, expected, tolerance = 1e-7)

    start(expected[[1]]) <- 2394474 - 10500
    end(expected[[1]]) <- 2420377 - 34053
    mcols(expected[[1]])$rate <- 2350
    mcols(expected[[1]])$min_of_max <- 0
    mcols(expected[[1]])$min_of_avg <- 1
    pw <- polymeraseWave(
        reads[[1]], reads[[2]], genes, approxDist, progress=FALSE,
        emissionDistAssumption="norm", TSmooth="3RS3R")
    expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("polymeraseWave normal exponental distribution assumption", {
    end(expected[[1]]) <- 2420377 - 34553
    mcols(expected[[1]])$rate <- 700
    mcols(expected[[1]])$min_of_max <- 0
    mcols(expected[[1]])$min_of_avg <- 1
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
