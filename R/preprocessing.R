#' Read Tn5 Insertions from paired-end bam file
#'
#' @param bam.file
#' @param genome
#' @param is.tn5Adjusted
#' @param mapq
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom GenomicAlignments readGAlignmentsList readGAlignmentPairs first last
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom rtracklayer GRangesForBSGenome
#' @importFrom dplyr %>% group_by summarise n tbl_df
#'
#' @export
#'
readInsertions <- function(bam.file, genome, is.tn5Adjusted = F, mapq = 30) {

  ### Obtain genome size for given genome
  genome.size <- GRangesForBSGenome(genome = genome)
  chrs <- as.character(seqnames(genome.size))
  chrs <- chrs[-grep(chrs, pattern = '(_|M)')]
  genome.size <- genome.size[as.character(seqnames(genome.size)) %in% chrs,]
  seqlevels(genome.size) <- chrs

  ### Read *.bam file
  flag <- scanBamFlag(isUnmappedQuery = F,
                      isDuplicate = F,
                      isNotPassingQualityControls = F,
                      isPaired = T,
                      isProperPair = T)
  flags <- ScanBamParam(flag = flag, which = genome.size, what = c('mapq'))
  reads <- readGAlignmentsList(bam.file, param = flags, use.names = T)

  ### Subset reads into first and second fragment AND calculate fragmentSize
  reads <- unlist(reads)
  plus.reads <- reads[strand(reads) == '+']
  minus.reads <- reads[strand(reads) == '-']

  # Filter reads and keep only matching pairs
  if(length(plus.reads) != length(minus.reads)) {
    print('BAM file was not filtered for proper PE reads!')
    plus.reads <- plus.reads[names(plus.reads) %in% names(minus.reads)]
    minus.reads <- minus.reads[names(minus.reads) %in% names(plus.reads)]
  }

  # Order read pairs
  i.match <- match(names(plus.reads), names(minus.reads))
  plus.reads <- plus.reads[i.match,]

  #subs <- strand(first(reads)mean_n) == '+'
  #plus.reads <- first(reads)[subs]#, last(reads)[!subs])
  #minus.reads <- last(reads)[subs]#, first(reads)[!subs])
  # Compute insertions
  frag.chrs <- seqnames(plus.reads)
  frag.start <- start(plus.reads)
  frag.end <- end(minus.reads)
  # Compute fragment lengths
  frag.len <- frag.end - frag.start

  ### Adjust for Tn5 bias:
  if(!is.tn5Adjusted) {
    # Plus strand insertions
    frag.start <- frag.start + 4
    # Minus strand insertions
    frag.end <- frag.end - 5
    # Fragment lengths
    frag.len <- frag.len - 9
  }

  ### Filter MAPQ
  #i.rmv <- which(mcols(lreads)$mapq < mapq & mcols(rreads)$mapq < mapq)
  #if(length(i.rmv) > 0) print('NOT MAPQ FILTERED BAM!')

  ### Prepare GRanges with insertions
  # Plus strand
  insertions.gr <- GRanges(seqnames = frag.chrs,
                           ranges = IRanges(frag.start, frag.start),
                           strand = '+',
                           fraglen = frag.len)
  # Minus strand
  minus.gr <- GRanges(seqnames = frag.chrs,
                      ranges = IRanges(frag.end, frag.end),
                      strand = '-',
                      fraglen = frag.len)

  insertions.gr <- c(insertions.gr, minus.gr)
  return(insertions.gr)
}

#' Read Tn5 Insertions from single-end bam file
#'
#' @param bam.file
#' @param genome
#' @param is.tn5Adjusted
#' @param mapq
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom rtracklayer GRangesForBSGenome
#' @importFrom dplyr %>% group_by summarise n tbl_df
#'
#' @export
readInsertionsFromSE <- function(bam.file, genome, is.tn5Adjusted = F, mapq = 30) {

  ### Obtain genome size for given genome
  genome.size <- GRangesForBSGenome(genome = genome)
  chrs <- as.character(seqnames(genome.size))
  chrs <- chrs[-grep(chrs, pattern = '(_|M)')]
  genome.size <- genome.size[as.character(seqnames(genome.size)) %in% chrs,]
  seqlevels(genome.size) <- chrs

  ### Read *.bam file
  flag <- scanBamFlag(isUnmappedQuery = F,
                      isDuplicate = F,
                      isNotPassingQualityControls = F)
  flags <- ScanBamParam(flag = flag, which = genome.size, what = c('mapq'))
  reads <- readGAlignments(bam.file, param = flags, use.names = T)

  ### Split reads into +/- strand
  is.plus <- which(as.character(strand(reads)) == '+')
  read.start <- start(reads)[is.plus]
  read.end <- end(reads)[-is.plus]

  ### Adjust for Tn5 bias:
  if(!is.tn5Adjusted) {
    # Plus strand insertions
    read.start <- read.start + 4
    # Minus strand insertions
    read.end <- read.end - 5
  }

  ### Prepare GRanges with insertions
  # Plus strand
  insertions.gr <- GRanges(seqnames = seqnames(reads)[is.plus],
                           ranges = IRanges(read.start, read.start),
                           strand = '+')
  # Minus strand
  minus.gr <- GRanges(seqnames = seqnames(reads)[-is.plus],
                      ranges = IRanges(read.end, read.end),
                      strand = '-')

  # Combine Insertions to final GRange
  insertions.gr <- sort(c(insertions.gr, minus.gr))
  return(insertions.gr)
}


#' Import broad peak file
#'
#' @param peak.file
#'
#' @return
#'
#' @importFrom rtracklayer import
#'
#' @export
#'
importBroadPeak <- function(peak.file) {
  add.cols <- c(signal = 'numeric', pvalue = 'numeric',
                qvalue = 'numeric')
  import(peak.file, format = 'BED', extraCols = add.cols)
}

#' Import narrow peak file
#'
#' @param peak.file
#'
#' @return
#'
#' @importFrom rtracklayer import
#'
#' @export
#'
importNarrowPeak <- function(peak.file) {
  add.cols <- c(singnal = 'numeric', pvalue = 'numeric',
                qvalue = 'numeric', peak = 'integer')
  import(peak.file, format = 'BED', extraCols = add.cols)
}


