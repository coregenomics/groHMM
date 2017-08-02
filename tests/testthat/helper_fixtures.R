library(groHMM)

tx <- GRanges(c("chr7:1000-30000"), strand="+")
annox <- GRanges(c("chr7:1000-11000",
                   "chr7:20000-30000"), strand="+")

## Fixtures (identical to example in polymeraseWave man page)
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
