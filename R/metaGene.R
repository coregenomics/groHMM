###########################################################################
##
##   Copyright 2013, 2014 Charles Danko and Minho Chae.
##
##   This program is part of the groHMM R package
##
##   groHMM is free software: you can redistribute it and/or modify it
##   under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful, but
##   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
##   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
##   for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##########################################################################

#' Returns a histogram of the number of reads in each section of a moving
#' window centered on a certain feature.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.
#' To change the number of processors, set the option 'mc.cores'.
#'
#' @param features A GRanges object representing a set of genomic coordinates.
#' The meta-plot will be centered on the transcription start site (TSS)
#' @param reads A GRanges object representing a set of mapped reads.
#' Instead of 'reads', 'plusCVG' and 'minusCVG' can be used  Default: NULL
#' @param plusCVG A RangesList object for reads with '+' strand.
#' @param minusCVG A RangesList object for reads with '-' strand.
#' @param size The size of the moving window.
#' @param up Distance upstream of each features to align and histogram.
#' Default: 10 kb.
#' @param down Distance downstream of each features to align and histogram.
#' If NULL, same as up. Default: NULL.
#' @param ... Extra argument passed to mclapply
#' @return Returns a integer-Rle representing the 'typical' signal
#' centered on a point of interest.
#' @author Charles G. Danko and Minho Chae
#' @examples
#' library(GenomicRanges)
#' features <- GRanges("chr7", IRanges(1000, 1000), strand="+")
#' reads <- GRanges("chr7", IRanges(start=c(1000:1004, 1100),
#'  width=rep(1, 6)), strand="+")
#' mg <- metaGene(features, reads, size=4, up=10)
metaGene <- function(features, reads=NULL, plusCVG=NULL, minusCVG=NULL,
    size=100L, up=10000L, down=NULL, ...) {
    seqlevels(features) <- seqlevelsInUse(features)
    ## Check 'reads'
    if (is.null(reads)) {
        if (is.null(plusCVG) || is.null(minusCVG))
            stop("Either 'reads' or 'plusCVG' and 'minusCVG' must be used")
    } else {
        seqlevels(reads) <- seqlevelsInUse(reads)
        plusCVG <- coverage(reads[strand(reads)=="+", ])
        minusCVG <- coverage(reads[strand(reads)=="-", ])
    }
    if (is.null(down)) down <- up

    featureList <- split(features, seqnames(features))

    H <- mclapply(
        seqlevels(features), metaGene_foreachChrom, featureList=featureList,
        plusCVG=plusCVG, minusCVG=minusCVG, size=size, up=up, down=down, ...)
    M <- sapply(seq_len(length(H)), function(x) as.integer(H[[x]]))

    return(Rle(apply(M, 1, sum)))
}


metaGene_foreachChrom <- function(chrom, featureList, plusCVG, minusCVG,
    size, up, down) {
    f <- featureList[[chrom]]

    pCVG <- plusCVG[[chrom]]
    mCVG <- minusCVG[[chrom]]
    offset <- floor(size/2L)

    pro <- promoters(f, upstream=up+offset, downstream=down + offset - 1L)

    M <- sapply(1:length(pro), function(x) {
        if (as.character(strand(pro)[x]) == "+")
            as.integer(runsum(
                pCVG[start(pro)[x]:end(pro)[x]], k=size))
        else
            as.integer(rev(runsum(
                mCVG[start(pro)[x]:end(pro)[x]], k=size)))
    })
    return(Rle(apply(M, 1, sum)))
}

#' Runs meta gene analysis for sense and anti-sense direction.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.
#' To change the number of processors, set the option 'mc.cores'.
#'
#' @param features GRanges A GRanges object representing a set of genomic
#' coordinates, i.e., set of genes.
#' @param reads GRanges of reads.
#' @param anchorType Either 'TSS' or 'TTS'.  Meta gene will be centered on the
#' transcription start site(TSS) or transcription termination site(TTS).
#' Default: TSS.
#' @param size Numeric.  The size of the moving window. Default: 100L
#' @param normCounts Numeric.  Normalization vector such as average reads.
#' Default: 1L
#' @param up Numeric. Distance upstream of each feature to align and histogram.
#' Default: 1 kb
#' @param down Numeric. Distance downstream of each feature to align and
#' histogram.  If NULL, down is same as up. Default: NULL
#' @param sampling Logical.  If TRUE, sub-sampling of meta gene is used.
#' Default: FALSE
#' @param nSampling Numeric. Number of sub-sampling.  Default: 1000L
#' @param samplingRatio Numeric. Ratio of sampling for features.  Default: 0.1
#' @param ... Extra argument passed to mclapply.
#' @return A list of integer-Rle for sense and anti-sense.
#' @author Minho Chae
#' @examples
#' library(GenomicRanges)
#' features <- GRanges("chr7", IRanges(start=1000:1001, width=rep(1,2)),
#'  strand=c("+", "-"))
#' reads <- GRanges("chr7", IRanges(start=c(1000:1003, 1100:1101),
#'  width=rep(1, 6)), strand=rep(c("+","-"), 3))
#' ## Not run:
#' # mg <- runMetaGene(features, reads, size=4, up=10)
runMetaGene <- function(features, reads, anchorType="TSS", size=100L,
    normCounts=1L, up=10000L, down=NULL, sampling=FALSE, nSampling=1000L,
    samplingRatio=0.1, ...) {
    # Check 'anchorType'
    if (!anchorType %in% c("TSS", "TTS")) {
        stop("'anchorType' must be either 'TSS' or 'TTS'")
    }

    if (anchorType == "TSS") {
        f <- resize(features, width=1L, fix="start")
    } else if (anchorType == "TTS") {
        f <- resize(features, width=1L, fix="end")
    }

    if (is.null(down)) down <- up

    fRev <- f
    strand(fRev) <- rev(strand(f))

    plusCVG <- coverage(reads[strand(reads)=="+", ])
    minusCVG <- coverage(reads[strand(reads)=="-", ])

    ## Sense direction
    ## nocov start
    if (sampling) {
        sense <- samplingMetaGene(
            features=f, plusCVG=plusCVG, minusCVG=minusCVG, size=size, up=up,
            down=down, nSampling=nSampling, samplingRatio=samplingRatio, ...)
        ## nocov end
    } else {
        sense <- metaGene(
            features=f, plusCVG=plusCVG, minusCVG=minusCVG, size=size, up=up,
            down=down, ...)
        sense <- sense/length(features)
    }

    ## Anti-sense direction
    ## nocov start
    if (sampling) {
        antisense <- samplingMetaGene(
            features=fRev, plusCVG=plusCVG, minusCVG=minusCVG, size=size, up=up,
            down=down, nSampling=nSampling, samplingRatio=samplingRatio, ...)
        ## nocov end
    } else {
        antisense <- metaGene(
            features=fRev, plusCVG=plusCVG, minusCVG=minusCVG, size=size, up=up,
            down=down, ...)
        antisense <- antisense/length(features)
    }

    sense <- sense * normCounts
    antisense <- antisense * normCounts
    return(list(sense=sense, antisense=antisense))
}


samplingMetaGene <- function(features, plusCVG, minusCVG, size=100L, up=10000L,
    down=NULL, nSampling=1000L, samplingRatio=0.1, ...) {
    ## nocov start
    samplingSize <- round(length(features)*samplingRatio)

    metaList <- mclapply(1:length(features), function(x) {
        metaGene(features=features[x, ], plusCVG=plusCVG, minusCVG=minusCVG,
        size=size, up=up, down=down)
    }
    , ...)

    allSamples <- mclapply(1:nSampling, function(x) {
        inx <- sample(1:length(features), size=samplingSize, replace=TRUE)
        onesample <- metaList[inx]
        mat <- sapply(onesample, function(x) as.integer(x))
        Rle(apply(mat, 1, sum))
    }
    , ...)

    M <- sapply(allSamples, function(x) as.integer(x))
    return(Rle(apply(M, 1, median)) / samplingSize)
    ## nocov end
}
