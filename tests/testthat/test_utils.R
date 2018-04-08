context("Utilities")

test_that("windowAnalysis returns list", {
    result <- windowAnalysis(reads = tx, windowSize = 100)
    expect_is(result, "list")
    expect_is(result[[1]], "Rle")
})

test_that("windowAnalysis returns correct Rle value", {
    seqlevels(tx) <- c("chr7", "chr8")
    tx_78 <- c(tx, GRanges(c("chr8:1000-30000"), strand="-"))
    expected <- list(chr7 = Rle(c(0, 1, 100), c(9, 1, 290)),
                     chr8 = Rle(c(0, 1, 100), c(9, 1, 290)))
    result <- windowAnalysis(reads = tx_78, windowSize = 100)
    expect_equal(expected, result)
})

test_that("windowAnalysis subsets by chrom", {
    seqlevels(tx) <- c("chr7", "chr8")
    tx_78 <- c(tx, GRanges(c("chr8:1000-30000"), strand="-"))
    expected <- list(chr7 = Rle(c(0, 1, 100), c(9, 1, 290)))
    result <- windowAnalysis(reads = tx, windowSize = 100, chrom="chr7")
    expect_equal(expected, result)
})

test_that("windowAnalysis returns empty list for empty GRanges", {
    tx <- GRanges()
    expected <- list()
    result <- windowAnalysis(reads = tx, windowSize = 100)
    expect_equal(expected, result)
})

test_that("windowAnalysis throws error for out of bounds windowSize", {
    expect_error(windowAnalysis(reads = tx, windowSize = 0), "range")
    expect_error(windowAnalysis(reads = tx, windowSize = end(tx) + 1), "range")
})

test_that("windowAnalysis throws error for out of bounds stepSize", {
    expect_error(windowAnalysis(reads = tx, stepSize = 0), "range")
})

test_that("windowAnalysis allows empty seqlevels", {
    ## Add empty seqlevels to reads.
    require(TxDb.Hsapiens.UCSC.hg19.knownGene) # nolint
    hg19_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene) # nolint
    idx_has <- which(names(hg19_seqinfo) %in% names(seqinfo(tx)))
    idx <- union(idx_has, 1:length(hg19_seqinfo))
    seqinfo(tx) <- hg19_seqinfo[names(hg19_seqinfo)[idx]]
    ## Set warnings as errors to catch any mclapply errors.
    warn <- getOption("warn")
    on.exit(options(warn = warn))
    options(warn = 2)
    expect_is(windowAnalysis(reads = tx, windowSize = 100, strand="-"), "list")
})

test_that("limitToXkb returns GRanges", {
    expect_s4_class(limitToXkb(annox), "GRanges")
    expect_warning(limitToXkb(GRanges()))
})

test_that("limitToXkb truncates based on size and applies offset", {
    offset <- formals(limitToXkb)$offset
    size <- formals(limitToXkb)$size
    size_ok <- size - offset - floor(size / 2)
    annox_ <- resize(annox, width=c(size_ok, size + 10))
    expected <- resize(shift(annox, offset),
                       width=c(size_ok, size) - offset)
    result <- limitToXkb(annox_)
    expect_equal(expected, result)
})
