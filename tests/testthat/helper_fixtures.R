library(GenomicRanges)
library(groHMM)

tx <- GRanges(c("chr7:1000-30000"), strand="+")
annox <- GRanges(
    c("chr7:1000-11000", "chr7:20000-30000"), strand="+",
    GENE_ID=c("itchy", "scratchy"))

## Fixtures (identical to example in polymeraseWave man page)
genes <- GRanges("chr7:2394474-2420377:+", SYMBOL="CYP2W1", ID="54905")
load("sr1.rda")
reads <- sr1
