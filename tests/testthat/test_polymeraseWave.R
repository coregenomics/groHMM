context("Polymerase wave detection")

## Fixtures
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

test_that("polymeraseWave returns expected data.frame", {
  pw <- polymeraseWave(reads[[1]], reads[[2]], genes, approxDist)
  expect_s3_class(pw, "data.frame")
  expect_equal(pw, expected, tolerance = 1e-7)
})

test_that("countMappableReadsInInterval validates NonMap", {
  expect_error(polymeraseWave(reads[[1]], reads[[2]], genes,
                        approxDist, NonMap=TRUE), "NonMap")
})
