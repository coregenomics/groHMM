###########################################################################
##
##   Copyright 2009, 2010, 2011, 2012, 2013 Charles Danko and Minho Chae.
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

#' limitToXkb truncates a set of genomic itnervals at a constant, maximum size.
#'
#' @param features A GRanges object representing a set of genomic coordinates. 
#' The meta-plot will be centered on the start position.
#' @param offset Starts the interval from this position relative to the start 
#' of each genomic features.
#' @param size Specifies the size of the window.
#' @return Returns GRanges object with new genomic coordiates.
#' @author Minho Chae and Charles G. Danko
#' @examples
#' tx <- GRanges("chr7", IRanges(1000, 30000), strand="+")
#' newTX <- limitToXkb(tx)
##  This function limits a genomic range to a samll region relative to the 
## transcription site.
limitToXkb <- function(features, offset=1000, size=13000) {
    w <- width(features)

    # 1. do nothing for w < offset 
    # 2. offset < w and w < size
    small  <- (offset < w) & (w < size)
    if (any(small)) { 
        features[small,] <- flank(features[small,], -1*(w[small]-offset), 
            start=FALSE)
    }

    # 2. w > size
    big  <- w > size 
    if (any(big)) {
        features[big,] <- resize(features[big,], width=size)

        bigPlus <- big & as.character(strand(features))=="+"
        if (any(bigPlus)) 
            start(features[bigPlus,]) <- start(features[bigPlus,]) + offset 

        bigMinus <- big & as.character(strand(features))=="-"
        if (any(bigMinus)) 
            end(features[bigMinus,]) <- end(features[bigMinus,]) - offset 
    }

    return(features)
}
 
#' readBed Returns a GenomicRanges object constrcuted from the specified bed 
#' file.
#'
#' Bed file format is assumed to be either four column: seqnames, start, end, 
#' strand columns; or six column: seqnames, start, end, name, score, and strand.
#' Three column format is also possible when there is no strand information.
#'
#' Any additional arguments availiable to read.table can be specified.
#'
#' @param file Path to the input file.
#' @param ... Extra argument passed to read.table
#' @return Returns GRanges object representing mapped reads.
#' @author Minho Chae and Charles G. Danko.
readBed <- function(file, ...) {
    df <- read.table(file, ...)
        if(NCOL(df) == 3) {
                colnames(df) <- c("seqnames", "start", "end")
                df <- cbind(df, strand=Rle("*", NROW(df)))
        }
        if(NCOL(df) == 4) colnames(df) <- c("seqnames", "start", "end", 
                                            "strand")
        if(NCOL(df) == 6) colnames(df) <- c("seqnames", "start", "end", 
                                            "name", "score", "strand")
        return( GRanges(seqnames = Rle(df$seqnames), ranges = 
                IRanges(df$start, df$end), strand = Rle(strand(df$strand))))
}
