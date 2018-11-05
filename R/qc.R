#' Function to compute the proportion of fragments within a given range
#' of fragment lengths.
#'
#' @param frag.lengths
#' @param min.length
#' @param max.length
#'
#' @return
#'
#' @importFrom dplyr data_frame %>% filter summarise n tbl_df
#'
computeFragmentProp <- function(frag.lengths, min.length, max.length) {
  # Compute total number of fragments
  n.frags <- sum(frag.lengths$n)
  y <- max(frag.lengths$freq)
  y <- y + y*0.05
  # Compute fragment proportion in segment
  frag.prop <- frag.lengths %>%
    filter(frag_length >= min.length & frag_length <= max.length) %>%
    summarise(x = min.length + (max.length - min.length)/2,
              y = y,
              prop_frags = round(sum(n)/n.frags*100, 2))
  return(frag.prop)
}

#' Plot fragment length distribution for ATAC-seq PE data.
#'
#' @param insertions.gr
#' @param min.length
#' @param max.length
#'
#' @return
#'
#' @import GenomicRanges ggplot2 cowplot
#' @importFrom dplyr data_frame %>% filter mutate summarise n tbl_df bind_rows
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#'
plotFragmentLengthDist <- function(insertions.gr, min.length = 0, max.length = 500) {
  # Filter for fragment starts and get fragment lengths
  i <- which(as.character(strand(insertions.gr)) == '+')
  frag.lengths <- insertions.gr$fraglen[i]
  # Compute fragment length distribution
  frag.lengths <- data_frame('frag_length' = frag.lengths) %>%
    group_by(frag_length) %>%
    summarise(n = n()) %>%
    mutate(freq = n/sum(n))
  # Compute proportion of nuc-free fragments and mono-nuc fragments
  nucfree.prop <- computeFragmentProp(frag.lengths, 0, 100) %>%
    mutate(prop_frags = sprintf('Nuc free\nfragments\n%s%%', prop_frags))
  mononuc.prop <- computeFragmentProp(frag.lengths, 150, 250) %>%
    mutate(prop_frags = sprintf('Mononuc\nfragments\n%s%%', prop_frags))
  nuc.prop <- rbind(nucfree.prop, mononuc.prop)
  # Plot Fragment length plot
  # Buenstro et.al 2013:
  # "Nucleosome positioning. To generate the nucleosome-position data track,
  # we chose to split reads into various bins. Reads below 100 bp were considered
  # nucleosome free, reads between 180 and 247 bp were considered to be mononucleosomes,
  # reads between 315 and 473 bp were considered to be dinucleosomes,
  # and reads between 558 and 615 bp were considered to be trinucleosomes"
  nuc.thres <- c(100, 150, 247)
  gg.fragLengthDist <- frag.lengths %>%
    filter(frag_length > min.length & frag_length < max.length) %>%
    ggplot(aes(x = frag_length, y = freq)) + geom_line(col = 'red') +
    scale_x_continuous(limits = c(0, max.length)) +
    geom_vline(xintercept = nuc.thres, col = 'black', linetype = 'dashed') +
    geom_label_repel(data = nuc.prop, aes(x = x, y = y, label = prop_frags)) +
    xlab('Fragment length (bp)') + ylab('Frequency') +
    theme_bw() + theme_cowplot()
  return(gg.fragLengthDist)
}

#' Compute fraction of insertions in peaks
#'
#' @param insertions.gr
#' @param peaks
#'
#' @return
#'
#' @import GenomicRanges
#'
#' @export
#'
computeFIP <- function(insertions.gr, peaks) {
  n.ov <- countOverlaps(insertions.gr, peaks)
  fip <- sum(n.ov > 0)/length(insertions.gr)
  fip <- round(fip*100, 2)
  return(fip)
}

#' Plots Tn5 insertion metaprofile around TSS and compute enrichment score
#'
#' @param insertions.gr
#' @param txdb
#' @param up
#' @param down
#' @param win.size
#'
#' @return
#'
#' @import GenomicRanges ggplot2 cowplot
#' @importFrom caTools runmean
#' @importFrom dplyr data_frame %>% filter mutate summarise
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
plotTssEnrichment <- function(insertions.gr, txdb, up = 2000, down = 2000, win.size = 50) {

  # Get promoter regions from TxDB object
  promoter.gr <- promoters(txdb, upstream = up, downstream = down)
  i.rmv <- grepl(seqnames(promoter.gr), pattern = '(_|M)')
  promoter.gr <- promoter.gr[!i.rmv,]

  # Compute promoter enrichment
  promoter.enrichment <- computeInsertionProfile(insertions.gr, promoter.gr,
                                                 strand.spec = T) %>%
    mutate(n = n/10^3) %>%
    mutate(smooth_n = runmean(n, win.size))

  # Compute TSS enrichment score
  tss.score <- promoter.enrichment %>%
    summarise(min_n = min(smooth_n), max_n = max(smooth_n)) %>%
    mutate(score = round(max_n/min_n, 2)) %>%
    mutate(label = sprintf('TSS enrichment\nscore = %s', score))

  # Plot promoter meta profile
  gg.promoter <- promoter.enrichment %>%
    ggplot(aes(x = rel_pos, y = n)) + geom_point() +
    geom_line(aes(y = smooth_n), col = 'red', size = 1) +
    geom_text_repel(data = tss.score, aes(x = -down/2, y = max_n, label = label)) +
    geom_vline(xintercept = 0, col = 'red', linetype = 'dashed') +
    geom_vline(xintercept = 250, col = 'black', linetype = 'dashed') +
    xlab('Position relative to TSS (bp)') +
    ylab('Number of insertions'~(x10^3)) + theme_cowplot()

  return(gg.promoter)
}

#' Computes TSS enrichment score
#'
#' @param insertions.gr
#' @param txdb
#' @param up
#' @param down
#' @param win.size
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom caTools runmean
#' @importFrom dplyr data_frame %>% filter mutate summarise
#'
#' @export
#'
computeTssEnrichmentScore <- function(insertions.gr, txdb, up = 2000, down = 2000, win.size = 50) {

  # Get promoter regions from TxDB object
  promoter.gr <- promoters(txdb, upstream = up, downstream = down)
  i.rmv <- grepl(seqnames(promoter.gr), pattern = '(_|M)')
  promoter.gr <- promoter.gr[!i.rmv,]

  # Compute promoter enrichment
  promoter.enrichment <- computeInsertionProfile(insertions.gr, promoter.gr,
                                                 strand.spec = T) %>%
    mutate(n = n/10^3) %>%
    mutate(smooth_n = runmean(n, win.size))

  # Compute TSS enrichment score
  tss.score <- promoter.enrichment %>%
    summarise(min_n = min(smooth_n), max_n = max(smooth_n)) %>%
    mutate(score = round(max_n/min_n, 2))

  return(tss.score)
}

