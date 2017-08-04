library(groHMM)

tx <- GenomicRanges::GRanges(c("chr7:1000-30000"), strand="+")
annox <- GenomicRanges::GRanges(c("chr7:1000-11000",
                                  "chr7:20000-30000"), strand="+")

## Fixtures (identical to example in polymeraseWave man page)
genes <- GenomicRanges::GRanges("chr7", IRanges::IRanges(2394474,2420377),
                                strand="+", SYMBOL="CYP2W1", ID="54905")
read_bams <- function(files) {
    lapply(files, function(x) {
        file_ <- system.file("extdata", x, package="groHMM", lib.loc=.Library)
        as(GenomicAlignments::readGAlignments(file_), "GRanges")
    })
}
reads <- GRangesList(read_bams(c("S0mR1.bam", "S40mR1.bam")))
