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

#' makeConsensusAnnotations Makes a consensus annotation
#'
#' Makes a non-overlapping consensus annotation.  Gene annotations are often
#' overlapping due to #' multiple isoforms for a gene.
#' In consensus annotation, isoforms are first reduced so that only
#' redundant intervals are used to represent a genomic interval for a gene,
#' i.e., a gene id.
#' Remaining unresolved annotations are further reduced by truncating 3'
#' end of annotations.
#'
#' @param ar GRanges of annotations to be collapsed.
#' @param minGap Minimum gap between overlapped annotations after truncated.
#' Default: 1L
#' @param minWidth Minimum width of consensus annotations. Default: 1000L
#' @param column Column by which to group transcripts.
#' @param BPPARAM Registered backend for BiocParallel.
#' Default: BiocParallel::bpparam()
#' @return Returns GRanges object of annotations.
#' @author Minho Chae
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' tx <- transcripts(
#'     TxDb.Hsapiens.UCSC.hg19.knownGene,
#'     columns=c("gene_id", "tx_id", "tx_name"),
#'     filter=list(tx_chrom="chr7"))
#' ## Workaround Travis-CI issue 7052
#' if (!is.na(Sys.getenv("TRAVIS", NA)))
#'     BiocParallel::register(BiocParallel::SerialParam())
#' ca <- makeConsensusAnnotations(tx)
makeConsensusAnnotations <- function(ar, minGap=1L, minWidth=1000L,
    column="gene_id", BPPARAM=bpparam()) {
    ## Check for gene ID column
    if (! any(colnames(mcols(ar)) %in% column))
        stop("Missing gene ID column")
    ## Drop rows missing gene ID values
    missing <- elementNROWS(mcols(ar)[, column]) == 0
    if (any(missing)) {
        ar <- ar[!missing, ]
        warning(
            sum(missing),
            " ranges do not have gene ID and they are dropped")
    }
    ## Drop multiple gene ID values
    many <- elementNROWS(mcols(ar)[, column]) > 1
    if (any(many)) {
        ar <- ar[!many, ]
        warning(
            sum(many),
            " ranges have multiple gene ID and they are dropped")
    }

    ar_list <- split(ar, unlist(mcols(ar)[, column]))
    singles <- unlist(ar_list[elementNROWS(ar_list) == 1])
    isoforms <- ar_list[elementNROWS(ar_list) > 1]

    message("Reduce isoforms(", length(isoforms), ") ... ", appendLF=FALSE)
    isoforms <- GRangesList(bplapply(isoforms, function(x) {
        ## For mixed strands or chrom, choose the longest
        if ( (length(seqlevelsInUse(x)) > 1) ||
            (length(unique(strand(x))) > 1)) {
            result <- x[which.max(width(x)), column]
        } else {
            dx <- disjoin(x)
            mcols(dx)[, column] <- mcols(x)[1, column]
            olcnt <- countOverlaps(dx, x)
            ## Use the disjoint ranges covered more than once
            multi <- dx[olcnt > 1]
            if (length(multi) == 0) {
                ## For non-overlapping isoforms, choose the longest
                result <- x[which.max(width(x)), column]
            } else if (length(multi) == 1) {
                result <- multi
            } else {
                reduced <- reduce(multi)
                if (length(reduced) == 1)
                    result <- reduced
                else (length(reduced) > 1)
                    result <- reduced[which.max(width(reduced)), ]

            }
            mcols(result)[, column] <- mcols(x)[1, column]
        }
        return(result)
    }
    , BPPARAM=BPPARAM))
    isoforms <- unlist(isoforms)
    message("OK")

    ## Check redundancy
    isoforms <- removeRedundant(isoforms)
    singles <- removeRedundant(singles)

    o <- findOverlaps(singles, isoforms, type="equal")
    if (length(o) != 0)
        singles <- singles[-queryHits(o), ]

    o <- findOverlaps(singles, isoforms, type="within")
    if (length(o) != 0)
        singles <- singles[-queryHits(o), ]

    o <- findOverlaps(isoforms, singles, type="within")
    if (length(o) != 0)
        isoforms <- isoforms[-queryHits(o), ]

    noiso <- sort(c(isoforms, singles[, column]))
    message("Truncate overlapped ranges ... ", appendLF=FALSE)
    ## with different gene_ids
    while (!isDisjoint(noiso)) {
        ol <- findOverlaps(noiso, drop.self=TRUE, drop.redundant=TRUE)
        ol_gr <- GRangesList(lapply(1:length(ol), function(x) {
            sort(c(noiso[queryHits(ol)[x]], noiso[subjectHits(ol)[x]]))
        }))

        ## Truncate 3' end
        ol_gr <- unlist(endoapply(ol_gr, function(x) {
            if (as.character(strand(x[1, ])) == "+") {
                end(x[1, ]) <- start(x[2, ]) - minGap
                ## first range's end is truncated
            } else {
                start(x[2, ]) <- end(x[1, ]) + minGap
                ## second range's end is truncated
            }
            x
        }))

        ## Remove any ranges with duplicated names since they already
        ## adjusted in the previous call
        ol_gr <- ol_gr[!duplicated(names(ol_gr)), ]

        noiso <- noiso[-unique(c(queryHits(ol), subjectHits(ol))), ]
        ## update noiso
        noiso <- c(noiso, ol_gr)
    }
    message("OK")

    noiso <- noiso[width(noiso) >= minWidth, ]
    return(sort(noiso))
}

removeRedundant <- function(annox) {
    o <- findOverlaps(annox, drop.self=TRUE, type="equal", drop.redundant=TRUE)
    if (length(o) != 0) annox <- annox[-subjectHits(o), ]

    o <- findOverlaps(annox, drop.self=TRUE, type="within", drop.redundant=TRUE)
    if (length(o) != 0) annox <- annox[-queryHits(o), ]

    return(annox)
}
