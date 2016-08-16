### R code from vignette source 'groHMM.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()
savewd = getwd()


###################################################
### code chunk number 2: install groHMM (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("groHMM")


###################################################
### code chunk number 3: install packages (eval = FALSE)
###################################################
## biocLite("GenomicFeatures")
## biocLite("org.Hs.eg.db")
## biocLite("edgeR")
## biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")


###################################################
### code chunk number 4: use groHMM
###################################################
library(groHMM)
options(mc.cores=getCores(4))


###################################################
### code chunk number 5: load data
###################################################
S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
        package="groHMM")), "GRanges")
S0mR2 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
        package="groHMM")), "GRanges") # Use R1 as R2
S40mR1 <- as(readGAlignments(system.file("extdata", "S40mR1.bam",
        package="groHMM")), "GRanges")
S40mR2 <- as(readGAlignments(system.file("extdata", "S40mR1.bam",
        package="groHMM")), "GRanges") # Use R1 as R2


###################################################
### code chunk number 6: combine replicates
###################################################
# Combine replicates
S0m <- c(S0mR1, S0mR2)
S40m <- c(S40mR1, S40mR2)


###################################################
### code chunk number 7: write wiggle files (eval = FALSE)
###################################################
## writeWiggle(reads=S0m, file="S0m_Plus.wig", fileType="wig", strand="+",
##         reverse=FALSE)
## writeWiggle(reads=S0m, file="S0m_Minus.wig", fileType="wig", strand="-",
##         reverse=TRUE)
## 
## # For BigWig file:
## # library(BSgenome.Hsapiens.UCSC.hg19)
## # si <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
## # writeWiggle(reads=S0m, file="S0m_Plus.wig", fileType="BigWig", strand="+",#
## #               reverse=FALSE, seqinfo=si)
## 
## # Normalized wiggle files
## expCounts <- mean(c(NROW(S0m), NROW(S40m)))
## writeWiggle(reads=S0m, file="S0m_Plus_Norm.wig", fileType="wig", strand="+",
##            normCounts=expCounts/NROW(S0m), reverse=FALSE)


###################################################
### code chunk number 8: transcript calling
###################################################
Sall <- sort(c(S0m, S40m))
# hmmResult <- detectTranscripts(Sall, LtProbB=-200, UTS=5,
#                   threshold=1)
# Load hmmResult from the saved previous run 
load(system.file("extdata", "Vignette.RData", package="groHMM"))
txHMM <- hmmResult$transcripts


###################################################
### code chunk number 9: called transcripts
###################################################
head(txHMM)


###################################################
### code chunk number 10: use ucsc genes
###################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
kgdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

library(GenomicFeatures)
# For refseq annotations:
# rgdb <- makeTxDbFromUCSC(genome="hg19", tablename="refGene")
# saveDb(hg19RGdb, file="hg19RefGene.sqlite")
# rgdb <- loadDb("hg19RefGene.sqlite")
kgChr7 <- transcripts(kgdb, columns=c("gene_id", "tx_id", "tx_name"),
                            filter=list(tx_chrom = "chr7"))
seqlevels(kgChr7) <- seqlevelsInUse(kgChr7)


###################################################
### code chunk number 11: make consensus annotations
###################################################
# Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(kgChr7, keytype="gene_id",
        mc.cores=getOption("mc.cores"))
library(org.Hs.eg.db)
map <- select(org.Hs.eg.db,
    keys=unlist(mcols(kgConsensus)$gene_id),
        columns=c("SYMBOL"), keytype=c("ENTREZID"))
mcols(kgConsensus)$symbol <- map$SYMBOL
mcols(kgConsensus)$type <- "gene"


###################################################
### code chunk number 12: evaluation of transcript calling
###################################################
e <- evaluateHMMInAnnotations(txHMM, kgConsensus)
e$eval


###################################################
### code chunk number 13: tuning HMM
###################################################
tune <- data.frame(LtProbB=c(rep(-100,3), rep(-200,3), rep(-300,3)),
           UTS=rep(c(5,10,15), 3))

Fp <- windowAnalysis(Sall, strand="+", windowSize=50)
Fm <- windowAnalysis(Sall, strand="-", windowSize=50)

# evals <- mclapply(seq_len(9), function(x) {
#        hmm <- detectTranscripts(Fp=Fp, Fm=Fm, LtProbB=tune$LtProbB[x],
#                  UTS=tune$UTS[x])
#        e <- evaluateHMMInAnnotations(hmm$transcripts, kgConsensus)
#        e$eval
#        }, mc.cores=getOption("mc.cores"),  mc.silent=TRUE)
tune <- cbind(tune, do.call(rbind, evals))  # evals from the previous run


###################################################
### code chunk number 14: check tuning
###################################################
tune
which.min(tune$total)
which.min(tune$errorRate)


###################################################
### code chunk number 15: expressed genes
###################################################
getExpressedAnnotations <- function(features, reads) {
     fLimit <- limitToXkb(features)
     count <- countOverlaps(fLimit, reads)
     features <- features[count!=0,]
     return(features[(quantile(width(features), .05) < width(features))
    & (width(features) < quantile(width(features), .95)),])}
conExpressed <- getExpressedAnnotations(features=kgConsensus,reads=Sall)


###################################################
### code chunk number 16: fig2plot
###################################################
td <- getTxDensity(txHMM, conExpressed, mc.cores=getOption("mc.cores"))
u <- par("usr")
lines(c(u[1], 0, 0, 1000, 1000, u[2]), c(0,0,u[4]-.04,u[4]-.04,0,0),
      col="red")
legend("topright", lty=c(1,1), col=c(2,1), c("ideal", "groHMM"))
text(c(-500,500), c(.05,.5), c("FivePrimeFP", "TP"))
td


###################################################
### code chunk number 17: fig2
###################################################
td <- getTxDensity(txHMM, conExpressed, mc.cores=getOption("mc.cores"))
u <- par("usr")
lines(c(u[1], 0, 0, 1000, 1000, u[2]), c(0,0,u[4]-.04,u[4]-.04,0,0),
      col="red")
legend("topright", lty=c(1,1), col=c(2,1), c("ideal", "groHMM"))
text(c(-500,500), c(.05,.5), c("FivePrimeFP", "TP"))
td


###################################################
### code chunk number 18: Non-mammalian (eval = FALSE)
###################################################
## # mysql --user=genome --host=genome-mysql.cse.ucsc.edu ce10 -e \
## #  "select chrom, txStart, txEnd, name, exonCount, strand, name2 from refGene \
## #   where chrom not like chrom!='chrM' and cdsStart != cdsEnd" | tail -n +1 > refgene.bed
## 
## # G <- read.table("refgene.bed", header=TRUE, stringsAsFactors=FALSE, sep="\t")
## # ce <- GRanges(G$chrom, IRanges(G$txStart, G$txEnd), strand=G$strand, \
## #               access=G$name, gene_id=G$name2)
## # ceConsensus <- makeConsensusAnnotations(ce, keytype="gene_id", \
## #                mc.cores=getOption("mc.cores"))


###################################################
### code chunk number 19: repairing called transcripts
###################################################
bPlus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="+")
bMinus <- breakTranscriptsOnGenes(txHMM, kgConsensus, strand="-")
txBroken <- c(bPlus, bMinus)
txFinal <- combineTranscripts(txBroken, kgConsensus)


###################################################
### code chunk number 20: plotTxDensity_final (eval = FALSE)
###################################################
## tdFinal <- getTxDensity(txFinal, conExpressed, mc.cores=getOption("mc.cores"))


###################################################
### code chunk number 21: DE for trascripts
###################################################
# For called transcripts
library(edgeR)
txLimit <- limitToXkb(txFinal)
ctS0mR1 <- countOverlaps(txLimit, S0mR1)
ctS0mR2 <- countOverlaps(txLimit, S0mR2)
ctS40mR1 <- countOverlaps(txLimit, S40mR1)
ctS40mR2 <- countOverlaps(txLimit, S40mR2)

pcounts <- as.matrix(data.frame(ctS0mR1, ctS0mR2, ctS40mR1, ctS40mR2))
group <- factor(c("S0m", "S0m", "S40m", "S40m"))
lib.size <- c(NROW(S0mR1), NROW(S0mR2), NROW(S40mR1), NROW(S40mR2))
d <- DGEList(counts=pcounts, lib.size=lib.size, group=group)
d <- estimateCommonDisp(d)
et <- exactTest(d)
de <- decideTestsDGE(et, p=0.001, adjust="fdr")
detags <- seq_len(NROW(d))[as.logical(de)]
# Number of transcripts regulated at 40m
cat("up: ",sum(de==1), " down: ", sum(de==-1), "\n")


###################################################
### code chunk number 22: fig3plot
###################################################
plotSmear(et, de.tags=detags)
# 2 fold up or down
abline(h = c(-1,1), col="blue")


###################################################
### code chunk number 23: fig3
###################################################
plotSmear(et, de.tags=detags)
# 2 fold up or down
abline(h = c(-1,1), col="blue")


###################################################
### code chunk number 24: DE for annotated genes
###################################################
# For ucsc knownGenes
kgChr7 <- transcripts(kgdb, columns=c("gene_id", "tx_id", "tx_name"),
                            filter=list(tx_chrom="chr7"))
map <- select(org.Hs.eg.db,
        keys=unique(unlist(mcols(kgChr7)$gene_id)),
        columns=c("SYMBOL"), keytype=c("ENTREZID"))
missing <- elementNROWS(mcols(kgChr7)[,"gene_id"]) == 0
kgChr7 <- kgChr7[!missing,]

inx <- match(unlist(mcols(kgChr7)$gene_id), map$ENTREZID)
mcols(kgChr7)$symbol <- map[inx,"SYMBOL"]

kgLimit <- limitToXkb(kgChr7)
ctS0mR1 <- countOverlaps(kgLimit, S0mR1)
ctS0mR2 <- countOverlaps(kgLimit, S0mR2)
ctS40mR1 <- countOverlaps(kgLimit, S40mR1)
ctS40mR2 <- countOverlaps(kgLimit, S40mR2)

counts <- as.matrix(data.frame(ctS0mR1, ctS0mR2, ctS40mR1, ctS40mR2))
group <- factor(c("S0m", "S0m", "S40m", "S40m"))
lib.size <- c(NROW(S0mR1), NROW(S0mR2), NROW(S40mR1), NROW(S40mR2))
d <- DGEList(counts=counts, lib.size=lib.size, group=group)

d <- estimateCommonDisp(d)
et <- exactTest(d)
de <- decideTestsDGE(et, p=0.001, adjust="fdr")
detags <- seq_len(NROW(d))[as.logical(de)]
symbols <- mcols(kgChr7)$symbol
# Number of unique genes regulated at 40m
cat("up: ", NROW(unique(symbols[de==1])), "\n")
cat("down: ", NROW(unique(symbols[de==-1])), "\n")
plotSmear(et, de.tags=detags)
abline(h = c(-1,1), col="blue")


###################################################
### code chunk number 25: metagene analysis
###################################################
upGenes <- kgChr7[de==1,]
expReads <- mean(c(NROW(S0m), NROW(S40m)))
# Metagene around TSS
mg0m <- runMetaGene(features=upGenes, reads=S0m, size=100,
                       normCounts=expReads/NROW(S0m), sampling=FALSE,
               mc.cores=getOption("mc.cores"))
mg40m <- runMetaGene(features=upGenes, reads=S40m, size=100,
                        normCounts=expReads/NROW(S40m), sampling=FALSE,
                mc.cores=getOption("mc.cores"))


###################################################
### code chunk number 26: plot metagene
###################################################
plotMetaGene <- function(POS=c(-10000:+9999), mg, MIN, MAX){
    plot(POS, mg$sense, col="red", type="h", xlim=c(-5000, 5000),
        ylim=c(floor(MIN),ceiling(MAX)), ylab="Read Density",
        xlab="Position (relative to TSS)")
     points(POS, (-1*rev(mg$antisense)), col="blue", type="h")
     abline(mean(mg$sense[5000:8000]), 0, lty="dotted")
}

MAX <- max(c(mg0m$sense, mg40m$sense))
MIN <- -1*max(c(mg0m$antisense, mg40m$antisense))
plotMetaGene(mg=mg0m, MIN=MIN, MAX=MAX)
plotMetaGene(mg=mg40m, MIN=MIN, MAX=MAX)


###################################################
### code chunk number 27: metagene figures
###################################################
png(filename="metagene0m.png")
plotMetaGene(mg=mg0m, MIN=MIN, MAX=MAX)
dev.off()
png(filename="metagene40m.png", width=480, height=480)
plotMetaGene(mg=mg40m, MIN=MIN, MAX=MAX)
dev.off()


###################################################
### code chunk number 28: save data for vignette (eval = FALSE)
###################################################
## save(hmmResult, evals, file="Vignette.RData")


###################################################
### code chunk number 29: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 30: groHMM.Rnw:598-599
###################################################
setwd(savewd)


