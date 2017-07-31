library(groHMM)
context("Polymerase wave detection")

## Fixtures (identical to example in man page)
genes <- GRanges("chr7", IRanges(2394474,2420377), strand="+",
                 SYMBOL="CYP2W1", ID="54905")
non_map <- GRanges("chr7", IRanges(2394474,2420377), strand="+")
read_bams <- function(files) {
    lapply(files, function(x) {
        as(readGAlignments(system.file("extdata", x, package="groHMM")),
           "GRanges")
        })
}
    
reads <- read_bams(c("S0mR1.bam", "S40mR1.bam"))
        
approxDist <- 20000

expected <- data.frame(
  StartWave = 9350,
  EndWave = 33750,
  Rate = 24400,
  minOfMax = NA,
  minOfAvg = NA,
  ID = "CYP2W1",
  ExternalID = factor(54905)
)

polyWave <- function(...) {
  capture.output(suppressMessages(suppressWarnings(pw <- polymeraseWave(...))))
  pw
}

test_that("polymeraseWave returns expected data.frame", {
  pw <- polyWave(reads[[1]], reads[[2]], genes, approxDist)
  expect_s3_class(pw, "data.frame")
  expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("countMappableReadsInInterval validates NonMap", {
  expect_error(polyWave(reads[[1]], reads[[2]], genes,
                        approxDist, NonMap=TRUE), "NonMap")
})
