#' Compute Tn5 sequence preference from a subset of insertions
#'
#' @param insertions.gr
#' @param bs.genome
#' @param ext
#' @param chrs
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom Biostrings getSeq vcountPattern consensusMatrix
#'
#'
#' @export
#'
computeTn5Bias <- function(insertions.gr, bs.genome, ext = 10, chrs = c('chr16', 'chr19')) {
  # Get index for insertions on the given chrs
  i <- which(as.character(strand(insertions.gr)) == '+' &
             as.character(seqnames(insertions.gr)) %in% chrs)
  # Get sequences for insertions +/- extension
  ins.seq <- getSeq(bs.genome, insertions.gr[i] + ext)
  # Remove sequences containing 'N'
  i.rmv <- which(vcountPattern('N', ins.seq) != 0)
  # Compute Tn5 bias consensus matrix
  tn5.bias <- consensusMatrix(ins.seq[-i.rmv,])[1:4,]
  return(tn5.bias)
}

#' Compute insertions meta-profile centered on a given set of regions
#'
#' @param insertions.gr
#' @param regions
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom dplyr data_frame %>% mutate group_by summarise n if_else ungroup
#'
#' @export
#'
computeInsertionProfile <- function(insertions.gr, regions, strand.spec = F) {
  # Overlap insertions and regions
  ov <- findOverlaps(insertions.gr, regions, ignore.strand = T)
  # Compute insertion meta-profil
  meta.profil <- data_frame('ins_pos' = start(insertions.gr)[queryHits(ov)],
                            'mid' = mid(ranges(regions))[subjectHits(ov)],
                            'mid_strand' = as.character(strand(regions))[subjectHits(ov)]) %>%
    mutate(rel_pos = as.numeric(ins_pos - mid))
  # If strand specific - flip positions
  if(strand.spec) {
    meta.profil <- meta.profil %>%
      mutate(rel_pos = if_else(mid_strand == '-', -1*rel_pos, rel_pos))
  }
  meta.profil <- meta.profil %>%
    group_by(rel_pos) %>%
    summarise(n = n()) %>%
    mutate(mean_n = n/length(regions),
           freq = n/sum(n),
           ipm = (n*10^6)/length(insertions.gr))
  return(meta.profil)
}


#' Plots insertion meta-profile around a motif given a insertion profile object
#'
#' @param ins.profile
#' @param motif.w
#' @param col
#'
#' @return
#'
#' @import ggplot2 cowplot
#' @importFrom dplyr %>%
#'
#' @export
#'
plotInsertionProfile <- function(ins.profile, motif.w, col = 'black') {
  gg.insProfile <- ins.profile %>%
    ggplot(aes(x = rel_pos, y = freq)) + geom_line(size = 1, col = col) +
      geom_vline(xintercept = c(-1,1)*motif.w, col = 'black',
                 linetype = 'dashed', size = 0.75) +
      xlab('Position relative to motif (bp)') + ylab('Insertion frequency')
  return(gg.insProfile)
}

#' Compute V-plot insertion matrix centered on a given set of regions
#'
#' @param insertions.gr
#' @param regions
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom dplyr data_frame %>% mutate group_by summarise n if_else ungroup
#'
#' @export
#'
#' @examples
computeVPlotMatrix <- function(insertions.gr, regions) {
  # Overlap insertions and regions
  ov <- findOverlaps(insertions.gr, regions, ignore.strand = T)
  # Compute insertion meta-profil
  vplot.data <- data_frame('ins_pos' = start(insertions.gr)[queryHits(ov)],
                            'fraglen' = insertions.gr$fraglen[queryHits(ov)],
                            'mid' = mid(ranges(regions))[subjectHits(ov)],
                            'mid_strand' = as.character(strand(regions))[subjectHits(ov)]) %>%
    mutate(rel_pos = as.numeric(ins_pos - mid)) %>%
        #   rel_pos = if_else(mid_strand == '-', -1*rel_pos, rel_pos)) %>%
    group_by(rel_pos, fraglen) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(freq = n/sum(n)) %>%
    select(-n)
  return(vplot.data)
}

#' Plots V-Plot for a before compute V-Plot data object
#'
#' @param vplot.data
#' @param smooth.w
#' @param motif.w
#'
#' @return
#'
#' @import ggplot2 cowplot
#' @importFrom dplyr %>% select tbl_df mutate
#' @importFrom tidyr spread gather
#' @importFrom caTools runmean
#'
#' @export
#'
plotVPlot <- function(vplot.data, motif.w, smooth.w = 5) {
  # Re-format vplot.data to a matrix
  v.matrix <- vplot.data %>%
    spread(rel_pos, freq) %>%
    select(-fraglen)

  # Smooth V-plot frequencies
  v.matrix <- t(apply(v.matrix, 1, function(r) runmean(r, smooth.w)))
  colnames(v.matrix) <- unique(vplot.data$rel_pos)

  # Plot V-plot
  gg.vplot <- v.matrix %>%
    tbl_df() %>%
    mutate(fraglen = unique(vplot.data$fraglen)) %>%
    gather(rel_pos, freq, 1:ncol(v.matrix)) %>%
    mutate(rel_pos = as.numeric(rel_pos)) %>%
    ggplot(aes(x = rel_pos, y = fraglen, fill = freq*10^4)) + geom_tile() +
    scale_fill_gradient(low = 'white', high = 'red', name = 'Frequency'~(x10^-4)) +
      geom_vline(xintercept = c(-1,1)*motif.w, col = 'black',
                 linetype = 'dashed', size = 0.75) +
    xlab('Position relative to motif (bp)') + ylab('Fragment length (bp)') +
    theme(legend.position = 'bottom')
  return(gg.vplot)
}


#' Computes CENTIPEDE insertions matrix around given motif sites
#'
#' @param insertions.gr
#' @param motifs
#' @param motif.ext
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom dplyr data_frame mutate if_else group_by summarise n bind_rows select
#' @importFrom tidyr spread
#'
#' @export
#'
computeCentipedeMatrix <- function(insertions.gr, motifs, motif.ext = 100) {

  # Overlap motifs and insertions
  ov <- findOverlaps(insertions.gr, motifs + motif.ext, ignore.strand = T)

  # Compute insertions per motif
  motif.ins <- data_frame('motif_id' = subjectHits(ov),
                          'motif_mid' = mid(ranges(motifs))[subjectHits(ov)],
                          'motif_strand' = as.character(strand(motifs))[subjectHits(ov)],
                          'ins_pos' = start(insertions.gr)[queryHits(ov)],
                          'strand' = as.character(strand(insertions.gr))[queryHits(ov)],
                          'fraglen' = insertions.gr$fraglen[queryHits(ov)]) %>%
    mutate('rel_pos' = as.numeric(ins_pos - motif_mid),
           'rel_pos' = if_else(motif_strand == '+', rel_pos, -1*rel_pos)) %>%
    group_by(motif_id, rel_pos, strand) %>%
    summarise(n = n())

  # Add motifs with zero coverage
  miss.ids <- which(!1:length(motifs)%in%unique(motif.ins$motif_id))
  if(length(miss.ids) > 0) {
    motif.ins <- data_frame(motif_id = miss.ids, rel_pos = 0,
                            strand = '+', n = 0) %>%
      bind_rows(motif.ins, .)
  }

  # Get important motif parameters
  motif.width <- width(motifs)[1]
  motif.hf <- floor(motif.width/2)
  region.size <- 2*motif.ext + motif.width
  region.hf <- motif.ext + motif.hf

  # Prepare data for CENTIPEDE - convert data.frame to matrix
  motif.matrix <- motif.ins %>%
    ungroup() %>%
    mutate(rel_pos = rel_pos + region.hf,
           rel_pos = if_else(strand == '+', rel_pos, rel_pos + region.size + 1)) %>%
    dplyr::select(-strand) %>%
    tidyr::spread(rel_pos, n, fill = 0) %>%
    ungroup() %>%
    dplyr::select(-motif_id) %>%
    as.matrix()
  return(motif.matrix)
}

#' Compute ROC statistics for TF footprinting vs ChIP-seq comparison
#'
#' @param prob
#' @param labels
#' @param name
#'
#' @return
#'
#' @import ROCR
#' @import dplyr
#'
#' @export
#'
computeRocStats <- function(prob, labels, name = '') {
  # Compute ROC curve
  pred <- prediction(prob, labels)
  perf <- performance(pred, measure = 'tpr', x.measure = 'fpr')

  # Compute AUC
  auc <- performance(pred, measure = 'auc')
  auc <- round(auc@y.values[[1]], 2)

  # Prepare ROC stats output
  roc.stats<- data_frame(fpr = unlist(perf@x.values),
                         tpr = unlist(perf@y.values),
                         label = sprintf('%s (AUC = %s)', name, auc))
  return(roc.stats)
}

#==============================================================================#
#                             FUTURE FUNCTIONS
#==============================================================================#
#### PLOT TF footprint for bound/unbound given post. prob. threshold
#post.prob <- centFit$PostPr[,1] %>%
#  tbl_df() %>%
#  mutate(motif_id = 1:length(motifs)) %>%
#  rename(post_prob = value)
#
## Filter motifs by insertion threshold
#cov.thres <- 6 # 6 insertions == 3 fragments
#i.rmv <- which(apply(motif.matrix, 1, sum) < cov.thres)
#bound.motifs <- motifs[-i.rmv,]
#bound.motifs$prob <- post.prob$post_prob[-i.rmv]
#bound.motifs <- bound.motifs[order(bound.motifs$prob, decreasing = T),]
#
#motif.hf <- width(motifs)[1]/2
#heat.dist <- 125 - motif.hf
#c <- c(-1, 1) * (heat.dist + motif.hf)
#cov.heatmap <- heatmaps::CoverageHeatmap(windows = bound.motifs + heat.dist,
#                                         track = insertions.gr,
#                                         label = '', coords = c)
#
#clusters <- as.numeric(bound.motifs$prob > postProb.thres)
#footprint.heatmap <- image(cov.heatmap)
##i.rmv <- which(apply(footprint.heatmap, 1, sum) == 0)
##footprint.heatmap <- footprint.heatmap[-i.rmv,]
#footprint.heatmap <- t(apply(footprint.heatmap, 1, function(r) (r - mean(r))/sd(r)))
#footprint.heatmap <- heatmaps::Heatmap(footprint.heatmap, label = '', coords = c)
#smoothed.heatmap <- smoothHeatmap(footprint.heatmap)#, sigma = c(2, 2))
#scale(smoothed.heatmap) <- c(-1.5, 1.5)
#
##pdf(sprintf('%s_cuttingHeatmap.pdf', pwm.ids$tf), width = 5, height = 6)
#plotHeatmapList(smoothed.heatmap, color = c('blue', 'white', 'red'),
#                partition=c(sum(clusters==1), sum(clusters==0)),
#                partition.lines=TRUE, refline = TRUE, legend = TRUE,
#                legend.pos = 'l')
##dev.off()
