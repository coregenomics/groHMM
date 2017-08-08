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

#' limitToXkb truncates a set of genomic intervals at a constant, maximum size.
#'
#' @param features A GRanges object representing a set of genomic coordinates.
#' @param offset Starts the interval from this position relative to the start
#' of each genomic features.
#' @param size Specifies the size of the window.
#' @return Returns GRanges object with new genomic coordinates.
#' @author Minho Chae and Charles G. Danko
#' @examples
#' library(GenomicRanges)
#' tx <- GRanges("chr7", IRanges(1000, 30000), strand="+")
#' newTX <- limitToXkb(tx)
## This function limits a genomic range to a small region relative to the
## transcription site.
limitToXkb <- function(features, offset=1000, size=13000) {
    features <- .normArgRanges(features)
    w <- width(features)

    ## 1. do nothing for w < offset
    ## 2. offset < w and w < size
    small  <- (offset < w) & (w < size)
    if (any(small)) {
        features[small, ] <-
            flank(
                features[small, ], -1 * (w[small] - offset), start=FALSE)
    }

    ## 2. w > size
    big  <- w > size
    if (any(big)) {
        features[big, ] <- resize(features[big, ], width=size)

        bigPlus <- big & as.character(strand(features))=="+"
        if (any(bigPlus))
            start(features[bigPlus, ]) <- start(features[bigPlus, ]) + offset

        bigMinus <- big & as.character(strand(features))=="-"
        if (any(bigMinus))
            end(features[bigMinus, ]) <- end(features[bigMinus, ]) - offset
    }

    return(features)
}
