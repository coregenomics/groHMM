tx <- GRanges(c("chr7:1000-30000"), strand="+")
annox <- GRanges(c("chr7:1000-11000",
                   "chr7:20000-30000"), strand="+")
txC <- annox
annoxC <- tx
numTxBr <- function(annox) length(ranges(breakTranscriptsOnGenes(tx, annox)))
numTxCb <- function(txC) length(ranges(combineTranscripts(txC, annoxC)))
