context("Polymerase wave detection")

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
  pw <- polymeraseWave(reads[[1]], reads[[2]], genes, approxDist, progress=FALSE)
  expect_s3_class(pw, "data.frame")
  expect_equal(pw, expected, tolerance = 1e-7)
})
