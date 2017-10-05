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

#' Given GRO-seq data, identifies the location of the polymerase 'wave' in
#' up- or down- regulated genes.
#'
#' The model is a three state hidden Markov model (HMM).  States represent:
#' (1) the 5' end of genes upstream of the transcription start site,
#' (2) up-regulated sequence, and (3) the 3' end of the gene through the
#' polyadenylation site.
#'
#' The model computes differences in read counts between the two conditions.
#' Differences are assumed fit a functional form which can be specified by
#' the user (using the emissionDistAssumption argument).
#' Currently supported functional forms include a normal distribution
#' (good for GRO-seq data prepared using the circular ligation protocol),
#' a gamma distribution (good for 'spiky' ligation based GRO-seq data),
#' and a long-tailed normal+exponential distribution was implemented, but
#' never deployed.
#'
#' Initial parameter estimates are based on initial assumptions of
#' transcription rates taken from the literature.  Subsequently all parameters
#' are fit using Baum-Welch expectation maximization.
#'
#' Reference: Danko CG, Hah N, Luo X, Martins AL, Core L, Lis JT, Siepel A,
#' Kraus WL. Signaling Pathways Differentially Affect RNA Polymerase II
#' Initiation, Pausing, and Elongation Rate in Cells. Mol Cell.
#' 2013 Mar 19. doi:pii: S1097-2765(13)00171-8. 10.1016/j.molcel.2013.02.015.
#'
#' @param reads1 Mapped reads in time point 1.
#' @param reads2 Mapped reads in time point 2.
#' @param genes A set of genes in which to search for the wave.
#' @param approxDist The approximate position of the wave.
#' Suggest using 2000 [bp/ min] * time [min], for mammalian data.
#' @param size The size of the moving window. Suggest using 50 for direct
#' ligation data, and 200 for circular ligation data.  Default: 50.
#' @param upstreamDist The amount of upstream sequence to include
#' Default: 10 kb.
#' @param TSmooth Optionally, outlying windows are set a maximum value
#' over the inter-quantile interval, specified by TSmooth.
#' Reasonable value: 20; Default: NA (for no smoothing).  Users are encouraged
#' to use this parameter ONLY in combination with the normal distribution
#' assumptions.
#' @param emissionDistAssumption Takes values "norm", "normExp", and "gamma".
#' Specifies the functional form of the 'emission' distribution for states
#' I and II (i.e. 5' of the gene, and inside of the wave).
#' In our experience, "gamma" works best for highly-variable 'spiky' data,
#' and "norm" works for smooth data.  As a general rule of thumb, "gamma"
#' is used for libraries made using the direct ligation method, and "norm"
#' for circular ligation data.  Default: "gamma".
#' @param filterWindowSize Method returns 'quality' information for
#' each gene to which a wave was fit.  Included in these metrics are several
#' that define a moving window.  The moving window size is specified by
#' filterWindowSize.  Default: 10 kb.
#' @param progress Whether to show progress bar.  Default: TRUE
#' @param BPPARAM Registered backend for BiocParallel.
#' Default: BiocParallel::bpparam()
#' @return Returns list of GRanges with Pol II wave positions and any
#' BiocParallel errors caught.
#' @author Charles G. Danko
#' @examples
#' library(GenomicAlignments)
#' library(BiocParallel)
#' genes <- GRanges("chr7", IRanges(2394474,2420377), strand="+",
#'                  SYMBOL="CYP2W1", ID="54905")
#' reads1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
#'                              package="groHMM")), "GRanges")
#' reads2 <- as(readGAlignments(system.file("extdata", "S40mR1.bam",
#'                              package="groHMM")), "GRanges")
#' approxDist <- 2000*10
#' ## Often, HMMs fails to converge or distributions fail to fit. Therefore
#' ## don't stop on error.  SerialParam is being used here for building, but
#' ## MulticoreParam or SnowParam are of course more desirable in practice.
#' bpparam <- SerialParam(stop.on.error = FALSE)
#' bpresult <- polymeraseWave(reads1, reads2, genes, approxDist, BPPARAM=bpparam)
#' ## Collect successful fits.
#' gr <- unlist(List(bpresult[bpok(bpresult)]))
#' gr
##  Given GRO-seq data, identifies the location of the polymerase wave in up-
##  or down-regulated genes.  This version is based on a full Baum-Welch EM
##  implementation.
##
##  This is a three state HMM -- initial state representing the intergenic
##  region 5' of a gene, the second representing the initially up-regulated
##  region, and the third representing the remaining sequence of a gene.
##
##  We assume that upstream region is intergenic, and thus its emission
##  distribution is assumed to be a constant, set based on a representative
##  intergenic region.  This is accommodated in my [1,*) framework by keeping
##  the variance constant, and scaling the mean for each gene.
##
## Test with GREB1:     chr2:11,591,693-11,700,363
## GREB1 <- data.frame(chr="chr2", start=11591693, end=11700363, str="+")
polymeraseWave <- function(reads1, reads2, genes, approxDist, size = 50,
    upstreamDist = 10000, TSmooth=NA, emissionDistAssumption = "gamma",
    filterWindowSize = 10000, progress=TRUE,
    BPPARAM = BiocParallel::bpparam()) {
    ## Parameter checks.
    if ( (filterWindowSize / size) %% 1 != 0)
        stop("filterWindowSize must be a multiple of size")

    ## Subset genes to those with reads within the search distances.
    search_unsafe <- suppressWarnings(resize(
        promoters(genes, upstream=upstreamDist),
        width=width(genes) + upstreamDist + approxDist))
    search <- trim(search_unsafe)
    has_reads <-
        search %over% reads1 |
        search %over% reads2
    if (sum(has_reads) < NROW(genes))
        message(
            "Dropping ", NROW(genes) - sum(has_reads), " of ", NROW(genes),
            " genes with no reads within or nearby")
    too_close_to_edge <- search_unsafe != search
    if (any(too_close_to_edge))
        message(
            "Dropping ", sum(too_close_to_edge), " of ", NROW(genes),
            " genes too close to chromosome edge")
    genes <- genes[has_reads & ! too_close_to_edge]
    if (NROW(genes) == 0) {
        message("No genes left to check")
        return(list(GRanges()))
    }

    Fp1 <- windowAnalysis(reads = reads1, strand = "+", windowSize = size)
    Fp2 <- windowAnalysis(reads = reads2, strand = "+", windowSize = size)
    Fm1 <- windowAnalysis(reads = reads1, strand = "-", windowSize = size)
    Fm2 <- windowAnalysis(reads = reads2, strand = "-", windowSize = size)

    ## Run the model separately on each gene.
    if (progress) {
        ## nocov start
        pb <-
            progress::progress_bar$
            new(
                show_after=0, clear=FALSE, total=NROW(genes),
                format=paste0(
                    " HMM [:bar] :percent :current/:total",
                    " elapsed: :elapsed eta: :eta"))
        pb$tick(0)
    } ## nocov end

    bptry(bplapply(seq_along(genes), function(i) {
        ## Define the gene in terms of the windowed size.
        chrom <- as.character(seqnames(genes[i]))
        if (all(strand(genes[i]) == "+")) {
            start <- floor( (start(genes[i]) - upstreamDist) / size)
            end   <- ceiling(end(genes[i]) / size)
            emis1  <- Fp1[[chrom]][c(start:end)]
            emis2  <- Fp2[[chrom]][c(start:end)]
        }
        else {
            start <- floor(start(genes[i]) / size)
            end   <- ceiling( (end(genes[i]) + upstreamDist) / size)
            emis1  <- rev(Fm1[[chrom]][c(start:end)])
            emis2  <- rev(Fm2[[chrom]][c(start:end)])
        }

        ## Scale to a minimum of 1 read at each position (for fitting Gamma).
        gene  <- as.numeric(emis1 - emis2)
        if (emissionDistAssumption == "gamma") {
            ## Leave centered on 0 for the norm_exp/norm emission functions
            gene  <- gene + (-1) * min(gene) + 1
            ## Must translate points if gamma distributed (gamma undefined <0).
        }

        if (is.double(TSmooth)) {
            ## Interprets it as a fold over the inter quantile interval to
            ## filter.
            medGene <- median(gene)
            iqrGene <- IQR(gene)
            gene[(medGene - gene) > (TSmooth * (iqrGene + 1))] <-
                medGene - (TSmooth * (iqrGene + 1))
            gene[(gene - medGene) > (TSmooth * (iqrGene + 1))] <-
                medGene + (TSmooth * (iqrGene + 1))
        } else if (!is.na(TSmooth)) {
            gene  <- smooth(gene, kind = TSmooth)
        }

        ## Make the initial guess +5kb --> approxDist.
        uTrans <- as.integer(ceiling( (upstreamDist - 5000) / size))
        iTrans <- as.integer(ceiling( (upstreamDist + approxDist) / size))

        ## Run Baum-Welch
        ##
        ## Set up initial parameter estimates.
        ##
        ## Fit transition and initial probabilities.
        tProb  <- as.list(data.frame(
            log(c( (1 - 1 / uTrans), 1 / uTrans, 0)),
            log(c(0, (1 - 1 / (iTrans - uTrans)), 1 / (iTrans - uTrans))),
            log(c(0, 0, 1))))  # Trans. prob.
        iProb  <- as.double(log(c(1, 0, 0))) # iProb.

        ## Fit initial distribution parameters for emission probabilities.
        parInt  <- Rnorm(gene[c(1:uTrans)])
        ## Check that the variance of the intergenic state is NOT 0.
        if (is.na(parInt$var) | parInt$var == 0) parInt$var <- 0.00001

        n_gene <- length(gene)
        if (emissionDistAssumption == "norm") {
            ePrDist <- c("norm", "norm", "norm")
            parPsi  <- Rnorm(gene[(uTrans + 1):iTrans])
            parBas  <- Rnorm(gene[(iTrans + 1):n_gene])
            ePrVars <- data.frame(
                c(parInt$mean, sqrt(parInt$var), -1, -1),
                c(parPsi$mean, sqrt(parPsi$var), -1, -1),
                c(parBas$mean, sqrt(parBas$var), -1, -1))
        }
        else if (emissionDistAssumption == "normExp") {
            ePrDist <- c("norm", "normexp", "normexp")
            parPsi  <- Rnorm.exp(gene[(uTrans + 1):iTrans], tol = 1e-4)
            parBas  <- Rnorm.exp(gene[(iTrans + 1):n_gene], tol = 1e-4)
            ePrVars <- data.frame(
                c(parInt$mean, sqrt(parInt$var), -1, -1), parPsi$parameters,
                parBas$parameters)
        }
        else if (emissionDistAssumption == "gamma") {
            ePrDist <- c("norm", "gamma", "gamma")
            parPsi  <- RgammaMLE(gene[(uTrans + 1):iTrans])
            parBas  <- RgammaMLE(gene[(iTrans + 1):n_gene])
            ePrVars <- data.frame(
                c(parInt$mean, sqrt(parInt$var), -1),
                c(parPsi$shape, parPsi$scale, -1),
                c(parBas$shape, parBas$scale, -1))
        }
        else {
            stop(
                "emissionDistAssumption should be set to: 'norm', ",
                "'normExp', or 'gamma'.")
        }
        ## Now run the HMM.
        g <- list()
        g[[1]] <- gene
        ans <- list()
        ans <- tryCatch(.Call(
            "RBaumWelchEM", as.integer(3), g, as.integer(1),
            ePrDist, ePrVars, tProb, iProb, 0.01, c(TRUE, TRUE, TRUE),
            c(TRUE, TRUE, TRUE), as.integer(10), FALSE, PACKAGE = "groHMM"),
            error = function(e) e)

        if (progress) pb$tick()         # nocov

        ##  Update emis...
        if (length(ans) < 3) {
            ## nocov start
            ## An error will have a length of 2 (is this guaranteed?!).
            message("Error caught on the C side")
            message(ans)
            ### Will be a previous iteration of ans...?! Generated a new 'ans'
            ans[[1]] <- NA
            ans[[2]] <- NA
            ans[[3]] <- c(rep(0, n_gene - 2), 1, 2)
            ans[[4]] <- NA
            ans[[5]] <- NA
        } ## nocov end

        ansVitervi <- ans[[3]][[1]]
        DTs <- max(which(ansVitervi == 0))
        DTe <- max(which(ansVitervi == 1))

        ## Calculate fit quality metrics using moving window.
        windowScale <- filterWindowSize / size
        MovMean  <- rep(0, n_gene)
        MovMax   <- rep(0, n_gene)
        for (k in 2:(n_gene - 2)) {
            MovMean[k] <- mean(gene[
                max(k - windowScale, 1):
                min(k + windowScale, n_gene)], na.rm = TRUE)
            MovMax[k]  <-  max(gene[
                max(k - windowScale, 1):
                min(k + windowScale, n_gene)], na.rm = TRUE)
        }
        if ( (DTs >= 1) & (DTe > 1) & (DTs < DTe) &
            (DTe < n_gene) & (DTs < n_gene)) {
            ## Iff converges to something useful.
            ANS <- (DTe - DTs) * size
            STRTwave <- DTs * size
            ENDwave <- DTe * size
            ## Calculates min/max and min/avg filters.
            medDns <- median(MovMax[
                max(which(ansVitervi == 1) + round(windowScale)):
                length(MovMax)
            ])
            minMax <- min(MovMax[c(
                min(which(ansVitervi == 1)):
                max(which(ansVitervi == 1)))])
            minWindLTMed <- as.numeric(medDns < minMax)
            ## True if min(wave) > med(wave.upstream)
            avgDns <- median(MovMean[
                max(which(ansVitervi == 1) + round(windowScale)):
                length(MovMean)
            ])
            minAvg <- min(MovMean[
                min(which(ansVitervi == 1)):
                max(which(ansVitervi == 1))
            ])
            minMeanWindLTMed <- as.numeric(avgDns < minAvg)

            gr <- genes[i]
            gr <- shift(gr, STRTwave)
            width(gr) <- ENDwave - STRTwave + 1
            mcols(gr)$rate <- ANS
            mcols(gr)$min_of_max <- minWindLTMed
            mcols(gr)$min_of_avg <- minMeanWindLTMed
            gr
        } else {
            stop("HMM failed to converge on anything useful.")
        }
    }
    , BPPARAM=BPPARAM))
}
