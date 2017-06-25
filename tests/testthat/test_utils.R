library(groHMM)
context("Utilities")

# Fixtures
reads <- GRanges(c("chr7:1000-30000"), strand="+")

test_that("windowAnalysis returns list", {
  result <- windowAnalysis(reads = reads, windowSize = 100)
  expect_is(result, "list")
  lapply(result, expect_is, "Rle")
})

test_that("windowAnalysis returns correct Rle value", {
  expected <- list(chr7 = Rle(c(0, 1, 100), c(9, 1, 290)))
  result <- windowAnalysis(reads = reads, windowSize = 100)
  expect_equal(expected, result)
})

test_that("windowAnalysis throws error with empty GRanges", {
  reads <- GRanges()
  expect_error(windowAnalysis(reads = reads, windowSize = 100), "cannot be empty")
})

test_that("windowAnalysis allows empty seqlevels", {
  # Add empty seqlevels to reads.
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  hg19_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)
  idx_has <- which(names(hg19_seqinfo) %in% names(seqinfo(reads)))
  idx <- union(idx_has, 1:length(hg19_seqinfo))
  seqinfo(reads) <- hg19_seqinfo[names(hg19_seqinfo)[idx]]
  # Set warnings as errors to catch any mclapply errors.
  on.exit(options(warn = getOption("warn")))
  options(warn = 2)
  expect_is(windowAnalysis(reads = reads, windowSize = 100, strand="-"), "list")
})
