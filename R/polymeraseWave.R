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
#' (2) upregulated sequence, and (3) the 3' end of the gene through the
#' polyadenylation site.
#'
#' The model computes differences in read counts between the two conditions.
#' Differences are assumed fit a functional form which can be specified by
#' the user (using the emissionDistAssumption argument).
#' Currently supported functional forms include a normal distribution
#' (good for GRO-seq data prepared using the circular ligation protocol),
#' a gamma distribution (good for 'spikey' ligation based GRO-seq data),
#' and a long-tailed normal+exponential distribution was implemented, but
#' never deployed.
#'
#' Initial parameter estimates are based on initial assumptions of
#' transcription rates taken from the literature.  Subsequently all parameters
#' are fit using Baum-Welch expetation maximization.
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
#' @param TSmooth Optimonally, outlying windows are set a maximum value
#' over the inter-quantile interval, specified by TSmooth.
#' Reasonable value: 20; Default: NA (for no smoothing).  Users are encouraged
#' to use this parameter ONLY in combination with the normal distribution
#' assumptions.
#' @param emissionDistAssumption Takes values "norm", "normExp", and "gamma".
#' Specifies the functional form of the 'emission' distribution for states
#' I and II (i.e. 5' of the gene, and inside of the wave).
#' In our experience, "gamma" works best for highly-variable 'spikey' data,
#' and "norm" works for smooth data.  As a general rule of thumb, "gamma"
#' is used for libraries made using the direct ligation method, and "norm"
#' for circular ligation data.  Default: "gamma".
#' @param finterWindowSize Method returns 'quality' information for
#' each gene to which a wave was fit.  Included in these metrics are several
#' that define a moving window.  The moving window size is specified by
#' filterWindowSize.  Default: 10 kb.
#' @param progress Whether to show progress bar.  Default: TRUE
#' @param ... Extra argument passed to mclapply
#' @return Returns either a data.frame with Pol II wave end positions,
#' or a List() structure with additional data, as specified by returnVal.
#' @author Charles G. Danko
#' @examples
#' genes <- GRanges("chr7", IRanges(2394474,2420377), strand="+",
#'                  SYMBOL="CYP2W1", ID="54905")
#' reads1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
#'                              package="groHMM")), "GRanges")
#' reads2 <- as(readGAlignments(system.file("extdata", "S40mR1.bam",
#'                              package="groHMM")), "GRanges")
#' approxDist <- 2000*10
#' # Not run:
#' # pw <- polymeraseWave(reads1, reads2, genes, approxDist)
##  Given GRO-seq data, identifies the location of the polymerase wave in up-
##  or down-regulated genes.  This version is based on a full Baum-Welch EM
##  implementation.
##
##  This is a three state HMM -- initial state representing the intergenic
##  region 5' of a gene, the second representing the initially upregulated
##  region, and the third representing the remaining sequence of a gene.
##
##  We assume that upstream region is intergenic, and thus its emmission
##  distriubtion is assumed to be a constant, set based on a representative
##  intergenic region.  This is accomidated in my [1,*) framework by keeping
##  the vairence constant, and scaling the mean for each gene.
##
## Test with GREB1:     chr2:11,591,693-11,700,363
## GREB1 <- data.frame(chr="chr2", start=11591693, end=11700363, str="+")
polymeraseWave <- function(reads1, reads2, genes, approxDist, size = 50,
    upstreamDist = 10000, TSmooth=NA, emissionDistAssumption = "gamma",
    finterWindowSize = 10000, progress = TRUE, ...) {

    genes <- as.data.frame(genes)
    genes <- genes[,c("seqnames", "start", "end", "strand", "SYMBOL", "ID")]

    Fp1 <- windowAnalysis(reads = reads1, strand = "+", windowSize = size)
    Fp2 <- windowAnalysis(reads = reads2, strand = "+", windowSize = size)
    Fm1 <- windowAnalysis(reads = reads1, strand = "-", windowSize = size)
    Fm2 <- windowAnalysis(reads = reads2, strand = "-", windowSize = size)

    ## Run the model separately on each gene.
    if (progress)
      pb <- progress::progress_bar$new(
        total = NROW(genes), format = "  HMM [:bar] :percent eta: :eta",
        clear = FALSE)
    rows <- mclapply(seq_along(NROW(genes)), function(i) {
      if (progress)
        pb$tick()
        ## Define the gene in terms of the windowed size.
        if (genes[i,4] == "+") {
            start <- floor((genes[i,2] - upstreamDist) / size)
            end   <- ceiling(genes[i,3] / size)
            emis1  <- 
                (as.numeric(Fp1[[ as.character(genes[i,1]) ]]))[c(start:end)]
            emis2  <- 
                (as.numeric(Fp2[[ as.character(genes[i,1]) ]]))[c(start:end)]
        }
        else {
            start <- floor(genes[i,2]/size)
            end   <- ceiling((genes[i,3] + upstreamDist) / size)
            emis1  <- 
                rev((as.integer(Fm1[[ as.character(genes[i,1]) ]]))
                    [c(start:end)])
            emis2  <- 
                rev((as.integer(Fm2[[ as.character(genes[i,1]) ]]))
                    [c(start:end)])
        }
    
        ## Scale to a minimum of 1 read at each position (for fitting Gamma). 
        gene  <- as.numeric(emis1 - emis2)
        if (emissionDistAssumption == "gamma") { 
            ## Leave centered on 0 for the norm_exp/norm emission functions
          gene  <- gene + (-1)*(min(gene)) + 1 
          ## Must translate points if gamma distributed (gamma undefined <0).
        }
        
        if (is.double(TSmooth)) { ## Interperts it as a fold over the inter 
                                  ## quantile interval to filter.
            message("TSmooth is.integer:", TSmooth)
            medGene <- median(gene)
            iqrGene <- IQR(gene)
            gene[(medGene - gene) > (TSmooth*(iqrGene + 1))] <- 
                medGene - (TSmooth*(iqrGene + 1))
            gene[(gene - medGene) > (TSmooth*(iqrGene + 1))] <- 
                medGene + (TSmooth*(iqrGene + 1))
        } else if (!is.na(TSmooth)) {
           gene  <- smooth(gene, kind = TSmooth)
        }

        ## Make the initial guess +5kb --> approxDist.
        uTrans <- as.integer(ceiling((upstreamDist - 5000)/size))
        iTrans <- as.integer(ceiling((upstreamDist + approxDist)/size))

        ## Run Baum-Welch
        MovMeanSpd <- finterWindowSize#10000
        Means <- rep(0,NROW(gene))
        MovMean  <- rep(0,NROW(gene))
        MovMax   <- rep(0,NROW(gene))
        dMovMean <- rep(0,NROW(gene))

        for (k in c(2:(NROW(gene) - 2))) {
            left  <- gene[c(1:k)]
            right <- gene[c((k + 1):NROW(gene))]
            Means[k]  <- mean(left, na.rm = TRUE) - mean(right, na.rm = TRUE)

            MovMean[k] <- mean(gene[
              max((k - (MovMeanSpd/size)),1):min((k + (MovMeanSpd/size)),
                                                 NROW(gene))], na.rm = TRUE)
            MovMax[k]  <-  max(gene[
              max((k - (MovMeanSpd/size)),1):min((k + (MovMeanSpd/size)),
                                                 NROW(gene))], na.rm = TRUE)
            dMovMean[k] <- MovMean[k - 1] - MovMean[k]
        }

        ########################################################################
        # Set up initial paremeter estimates.

        ## Fit transition and initial probabilities.
        tProb  <- as.list(data.frame(
            log(c((1 - (1/uTrans)),(1/uTrans),0)),
            log(c(0,(1 - (1/(iTrans - uTrans))),(1/(iTrans - uTrans)))), 
            log(c(0, 0, 1))))  # Trans. prob.
        iProb  <- as.double(log(c(1, 0, 0))) # iProb.

        ## Fit initial distribution paremeters for emission probabilities.
        parInt  <- Rnorm(gene[c(1:uTrans)])
        if (is.na(parInt$var) | parInt$var == 0) parInt$var = 0.00001 
        ## Check that the varience of the intergenic state is NOT 0.

        if (emissionDistAssumption == "norm") {
            ePrDist <- c("norm", "norm", "norm")
            parPsi  <- Rnorm(gene[c((uTrans + 1):iTrans)])
            parBas  <- Rnorm(gene[c((iTrans + 1):NROW(gene))])
            ePrVars <- data.frame(c(parInt$mean, 
                sqrt(parInt$var), -1, -1), 
                c(parPsi$mean, sqrt(parPsi$var), -1, -1), 
                c(parBas$mean, sqrt(parBas$var), -1, -1))
        }
        else if (emissionDistAssumption == "normExp") {
            ePrDist <- c("norm", "normexp", "normexp")
            parPsi  <- Rnorm.exp(gene[c((uTrans + 1):iTrans)], tol = 1e-4) #
            parBas  <- Rnorm.exp(gene[c((iTrans + 1):NROW(gene))], tol = 1e-4) #
            ePrVars <- data.frame(c(parInt$mean, sqrt(parInt$var), -1, -1),
                    parPsi$parameters, parBas$parameters)
        }
        else if (emissionDistAssumption == "gamma") {
            ePrDist <- c("norm", "gamma", "gamma")
            parPsi  <- RgammaMLE(gene[c((uTrans + 1):iTrans)])
            parBas  <- RgammaMLE(gene[c((iTrans + 1):NROW(gene))])
            ePrVars <- data.frame(c(parInt$mean, sqrt(parInt$var), -1), 
                c(parPsi$shape, parPsi$scale, -1), 
                c(parBas$shape, parBas$scale, -1))
        }
        else {
          stop("emissionDistAssumption should be set to: 'norm', ",
               "'normExp', or 'gamma'.")
        }
        ## Now run the HMM.
        g <- list()
        g[[1]] <- gene
        ans <- list()
        ans <- tryCatch(.Call("RBaumWelchEM", as.integer(3), g, as.integer(1),
                ePrDist, ePrVars, tProb, iProb, 0.01, c(TRUE,TRUE,TRUE),
                c(TRUE, TRUE, TRUE), as.integer(10), FALSE, PACKAGE = "groHMM"),
                error = function(e) e)
                                        ##  Update emis...
        if (NROW(ans) < 3) {
            ## An error will have a length of 2 (is this guaranteed?!).  
            print("ERROR CAUGHT ON THE C SIDE")
            print(ans) 
            ### Will be a previous iteration of ans...?! Generated a new 'ans'
            ans[[1]] <- NA
            ans[[2]] <- NA
            ans[[3]] <- c(rep(0, NROW(gene) - 2), 1, 2)
            ans[[4]] <- NA
            ans[[5]] <- NA
        }
        
        ansVitervi <- ans[[3]][[1]]
        DTs <- max(which(ansVitervi == 0))
        DTe <- max(which(ansVitervi == 1))

        ### Calculate quality filters... Wrap up.
        if ((DTs >= 1) & (DTe > 1) & (DTs < DTe) &
            (DTe < NROW(gene)) & (DTs < NROW(gene))) { 
            ## iff converges to something useful.
            ANS <- (DTe - DTs)*size
            STRTwave <- DTs*size
            ENDwave <- DTe*size
            ## Calculates min/max and min/avg filters.
            medDns <- median(MovMax[c(max((which(ansVitervi == 1)) +
                    round(MovMeanSpd/size)):NROW(MovMax))])
            minMax <- min(MovMax[c(min(which(ansVitervi == 1)):max(which(ansVitervi == 1)))])
            minWindLTMed <- (medDns < minMax)
            ## True if min(wave) > med(wave.upstream)
            avgDns <- median(MovMean[c(max((which(ansVitervi == 1)) +
                    round(MovMeanSpd/size)):NROW(MovMean))])
            minAvg <- min(MovMean[c(min(which(ansVitervi == 1)):max(which(ansVitervi == 1)))])
            minMeanWindLTMed <- (avgDns < minAvg)
        }
        data.frame(
          StartWave = STRTwave,
          EndWave = ENDwave,
          Rate = ANS,
          minOfMax = minWindLTMed,
          minOfAvg = minMeanWindLTMed)
    }, ...)  # Extra args to mclapply.
    
    cbind.data.frame(
      rbind.data.frame(rows),
      ID = genes[,5],
      ExternalID = genes[,6])
}
