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

#' Function identifies expressed features using the methods introduced in 
#' Core, Waterfall, Lis; Science, Dec. 2008.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors use the argument 'mc.cores'.
#'
#' @param features A GRanges object representing a set of genomic coordinates.
#' The meta-plot will be centered on the start position.  
#' @param reads A GRanges object representing a set of mapped reads.
#' @param Lambda Measurement of assay noise.  Default: 0.04 reads/ kb in a 
#' library of 10,751,533 mapped reads. (background computed in Core, 
#' Waterfall, Lis. (2008) Science.).
#' @param ... Extra argument passed to mclapply
#' @return Returns a data.frame representing the expression p.values for 
#' features of interest.
#' @author Charles G. Danko 
expressedGenes <- function(features, reads, Lambda=NULL, ...) {
    ## Order -- Make sure, b/c this is one of our main assumptions.  Otherwise
    ## violated for DBTSS.
    reads <- .normArgRanges(reads)
    reads <- reads[order(as.character(seqnames(reads)), start(reads)),] 
    C <- sort(unique(as.character(seqnames(features))))
    if(is.null(Lambda))
        ## NROW(reads) / genomeSize
        Lambda <- 0.04 * NROW(reads) / 10751533 / 1000
    
    ## Run parallel version.
    mcp <- mclapply(
        seq_along(C), expressedGenes_foreachChrom, C=C, features=features,
        reads=reads, Lambda=Lambda, ...)

    ## Unlist... 
    ANSpvalue <- rep(0,NROW(features))
    ANScounts <- rep(0,NROW(features))
    ANSgsize  <- rep(0,NROW(features))
    for(i in seq_along(C)) {
        indxF   <- which(as.character(seqnames(features)) == C[i])
        indxPrb   <- which(as.character(seqnames(reads)) == C[i])
        if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
            ANSpvalue[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["ANSpvalue"]]
            ANScounts[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["ANScounts"]]
            ANSgsize[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["ANSgsize"]]
        }
    }

    return(
        cbind.data.frame(
            data.frame(
                pval=ANSpvalue, readCounts=ANScounts, size=ANSgsize),
            mcols(features)))
}

expressedGenes_foreachChrom <- function(i, C, features, reads, Lambda) {
    ## Which KG?  prb?
    indxF   <- which(as.character(seqnames(features)) == C[i])
    indxPrb <- which(as.character(seqnames(reads)) == C[i])

    if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
        ## Order -- Make sure, b/c this is one of our main assumptions.  
        ## Otherwise violated for DBTSS.
        Ford <- order(start(features[indxF,])) 

        ## Type coersions.
        FeatureStart    <- start(features[indxF,][Ford])
        FeatureEnd  <- end(features[indxF,][Ford])
        FeatureStr  <- as.character(strand(features[indxF,][Ford]))
        PROBEStart  <- start(reads[indxPrb,])
        PROBEEnd    <- end(reads[indxPrb,])
        PROBEStr    <- as.character(strand(reads[indxPrb,]))

        ## Set dimensions.
        dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
        dim(FeatureEnd)     <- c(NROW(FeatureEnd),   NCOL(FeatureEnd))
        dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
        dim(PROBEStart)     <- c(NROW(PROBEStart),   NCOL(PROBEStart))
        dim(PROBEEnd)       <- c(NROW(PROBEEnd),     NCOL(PROBEEnd))
        dim(PROBEStr)       <- c(NROW(PROBEStr),     NCOL(PROBEStr))

        NUMReads <- .Call(
            "CountReadsInFeatures", FeatureStart, 
            FeatureEnd, FeatureStr, PROBEStart, PROBEEnd, 
            PROBEStr, PACKAGE="groHMM")

        ## Calculate Poisson prob. of each.
        ANSgsize_c <- FeatureEnd - FeatureStart
        ANSpvalue_c <- ppois(
            NUMReads, Lambda * ANSgsize_c, lower.tail=FALSE)
        ANScounts_c <- NUMReads

            
        return(
            list(
                ANSpvalue=ANSpvalue_c, ANScounts=ANScounts_c,
                ANSgsize=ANSgsize_c, ord=Ford))
    }
    return(integer(0))
}
