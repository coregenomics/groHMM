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

#' windowAnalysis Returns a vector of integers representing the counts of 
#' reads in a moving window.
#'
#' @param reads GenomicRanges object representing the position of reads 
#' mapping in the genome.
#' @param strand Takes values of "+", "-", or "*".  "*" denotes collapsing 
#' reads on both strands.  Default: "*".
#' @param windowSize Size of the moving window. Either windowSize or 
#' stepSize must be specified.
#' @param stepSize The number of bp moved with each step.
#' @param chrom Chromosome for which to return data.  
#' Default: returns all avaliable data.
#' @return Returns a list object, each element of which represents a 
#' chromosome.
#' @author Charles G. Danko and Minho Chae
#' @examples
#' library(GenomicAlignments)
#' S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
#'      package="groHMM")), "GRanges")
#' ## Not run:
#' # Fp <- windowAnalysis(S0mR1, strand="+", windowSize=50)
windowAnalysis <- function(reads, strand="*", windowSize=stepSize, 
    stepSize=windowSize, chrom=NULL) {
    reads <- .normArgRanges(reads, warnOnEmpty=FALSE)

    if (length(reads) == 0)
        return(list())

    if (!(windowSize > 0 & (windowSize <= max(end(reads)))))
        stop("'windowSize' is out of range!")

    if (! stepSize > 0)
        stop("'stepSize' is out of range!")

    if (!is.null(chrom))  
        reads <- reads[seqnames(reads) == chrom,]

    seqlevels(reads) <- seqlevelsInUse(reads)
    readsList <- split(reads, seqnames(reads))

    # Change reads' strand 
    readsList <- endoapply(readsList, function(x) {
            if (strand == "*") 
                strand(x) <- "*"
            else
                x <- x[strand(x) == strand,]
            x 
    })

    lapply(readsList, function(x) {
        cov <- coverage(x)[[1]]
        to <- (length(cov) %/% windowSize)*windowSize
        starts <- seq(1, to, stepSize)
        vi <- Views(cov, start=starts, width=windowSize)
        Rle(viewSums(vi))
    })
}
