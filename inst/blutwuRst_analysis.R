#' ---
#' title: 'Example analysis of ATAC-seq data using blutwuRst'
#' author: 'Sebastian Steinhauser'
#' date: '25/07/2017'
#' output:
#'   html_document:
#'     number_sections: yes
#'     toc: true
#'     toc_float: true
#'     fig_caption: yes
#'     code_folding: hide
#'   pdf_document:
#'     number_sections: yes
#'     toc: true
#'     fig_caption: yes
#' ---

#/*==========================================================================#*/
#' # Libraries, paths and parameters
#+ chunk_preparation, cache=F, error=F, warning=F, message=F, echo=T
#/*==========================================================================#*/
library(blutwuRst)
library(motifmatchr)
library(CENTIPEDE)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(dplyr)
library(tidyr)
library(readr)

library(ggplot2)
library(cowplot)

# Define path to bam/peak file
package.path <- '/Users/steinhs/scripts/R/packages/blutwuRst/'
bam.file <- file.path(package.path, 'data/D0_rmbqr_rmdup_sorted_chr19.bam')
peak.file <- file.path(package.path, 'data/D0_peaks_chr19.narrowPeak')

getwd()

# Define genome and transcriptome
genome <- 'mm10'
genome.gr <- rtracklayer::GRangesForBSGenome(genome)
genome.gr <- genome.gr[-grep(as.character(seqnames(genome.gr)), pattern ='(_|M)')]

bs.genome <- BSgenome.Mmusculus.UCSC.mm10

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#/*==========================================================================#*/
#' # Load insertions and peaks
#+ chunk_loadData, cache=F, error=F, warning=F, message=F, echo=T
#/*==========================================================================#*/
# Read insertions, if insertion.Rds not exists
ins.rds <- file.path(package.path, 'data/insertions.Rds')
if (!file.exists(ins.rds)) {
  insertions.gr <- readInsertions(bam.file, genome)
  write_rds(insertions.gr, ins.rds)
} else {
  insertions.gr <- read_rds(ins.rds)
}

# Load ATAC-seq peaks
atac.peaks <- importNarrowPeak(peak.file)

#/*==========================================================================#*/
#' # Quality Control
#+ chunk_qc, cache=F, error=F, warning=F, message=F, echo=T
#/*==========================================================================#*/
# Plot fragment length distribution
plotFragmentLengthDist(insertions.gr)

# Compute fraction of insertions in peaks
computeFIP(insertions.gr, atac.peaks)

# Compute TSS enrichment score & plot TSS enrichment
plotTssEnrichment(insertions.gr, txdb)

computeTssEnrichmentScore(insertions.gr, txdb)

#/*==========================================================================#*/
#' # TF footprinting
#+ chunk_tfFootprinting, cache=F, error=F, warning=F, message=F, echo=T
#/*==========================================================================#*/
library(JASPAR2016)
library(TFBSTools)

# Query CTCF motif
opts <- list(name = 'CTCF')
pwm.list <- getMatrixSet(JASPAR2016, opts)

# Find CTCF motifs on chr19s
motifs <- matchMotifs(pwms = pwm.list,
                      subject = genome.gr[seqnames(genome.gr) == 'chr19',],
                      genome = bs.genome,
                      out = 'position')

# Convert to regions +/-250 around motif
motif.regions <- motifs[[1]] + 250

# Compute insertion matrix
ins.profile <- computeInsertionProfile(insertions.gr[insertions.gr$fraglen <= 100,], motif.regions, strand.spec = T)

# Plot insertions around motif
ins.profile %>%
  ggplot(aes(x = rel_pos, y = freq)) + geom_line() +
  ylab('Tn5 insertion frequncy') + xlab('Position relativ to motif')

### RUN centipede
# Prepare annotation information for CENTIPEDE
motifIns.anno <- data_frame(x = 1, score = motifs[[1]]$score)

# Prepare inseration matrix using only nuc free insertions; +/- 100bp
ins.matrix <- computeCentipedeMatrix(insertions.gr[insertions.gr$fraglen <= 100], motifs[[1]])

# Fit CENTIPEDE model
cent.fit <- fitCentipede(Xlist = list('all' = ins.matrix),
                         Y = as.matrix(motifIns.anno))

plotProfile(cent.fit$LambdaParList[[1]], Mlen = width(motifs[[1]])[1])

