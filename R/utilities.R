#' Plot ggplot to multiple files defined by figure suffix
#'
#' @param gg.plot
#' @param file.name
#' @param figure.suffix
#' @param base_height
#' @param base_width
#'
#' @return
#'
#' @import ggplot2
#' @import cowplot
#'
#' @export
#'
savePlots <- function(gg.plot, file.name, figure.suffix = c('pdf', 'eps'), base_height = 4, base_width = NULL) {
  # Prepare outfile names
  out.files <- paste(file.name, figure.suffix, sep = '.')
  # Save plots to files
  lapply(out.files, function(out.file) {
    save_plot(plot = gg.plot, filename = out.file, base_height = base_height, base_width = base_width)
  })
}

#' Get summit position for GenomicRanges
#'
#' @param regions.gr
#'
#' @return
#'
#' @import GenomicRanges
#'
#' @export
#'
getSummit <- function(regions.gr) {
  # Get peak summits
  summits <- mid(ranges(regions.gr))
  # Prepare GRange from summits
  summits <- GRanges(seqnames = seqnames(regions.gr),
                     ranges = IRanges(summits, summits))
  return(summits)
}

#' Compute standard error of mean
#'
#' @param x
#'
#' @return
#'
#'
stdError <- function(x) {
  std.error <- sd(x)/sqrt(length(x))
  return(std.error)
}

#' Better version of listFiles
#'
#' @param path
#' @param is.file
#' @param pattern
#' @param maxdepth
#'
#' @return
#'
#' @export
#'
listFiles <- function(path, is.file = T, pattern = NULL, maxdepth = NULL) {
  # Paste cmd line from vars
  cmd <- path #sprintf('find %s', path)
  if(!is.null(maxdepth)) cmd <- paste(cmd, sprintf('-maxdepth %s', maxdepth))
  if(is.file) cmd <- paste(cmd, '-type f')
  if(!is.null(pattern)) cmd <- paste(cmd, sprintf('| grep %s', pattern))
  #files <- system(cmd, intern = T)
  files <- system2(command = 'find', args = cmd, stdout = T)
  return(files)
}

#' Save a GRange object as bed file
#'
#' @param grange
#' @param bed.file
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom dplyr %>% tbl_df mutate
#' @importFrom readr write_tsv
#'
#' @export
#'
grange2Bed <- function(grange, bed.file) {
  ranges <- as.data.frame(grange)[,1:3] %>%
    tbl_df() %>%
    mutate(name = 1:length(grange))
  write_tsv(x = ranges, path = bed.file, col_names = F)
}


#' Apply bigWigAverageOverBed to a GRange object
#'
#' @param grange
#' @param bw.file
#' @param is.minMax
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom readr read_tsv
#' @importFrom dplyr %>% arrange
#'
#' @export
#'
bwAverage4GRange <- function(grange, bw.file, is.minMax = T) {

  # Create tmp dir
  tmp.path <- file.path(tempdir(), 'bwAverage4GRange')
  if(!dir.exists(tmp.path)) dir.create(tmp.path, recursive = T)

  # Write GRange to bed file
  bed.file <- file.path(tmp.path, 'ranges.bed')
  grange2Bed(grange, bed.file)

  # Define output paht
  out.file <- file.path(tmp.path, 'out.txt')

  # Find bigWigAverageOverBed on Path variable and create cmd command
  cmd <- Sys.which(c('bigWigAverageOverBed'))
  if(length(cmd) == 0) warning('bigWigAverageOverBed is not on $PATH')
  cmd <- paste(cmd, bw.file, bed.file, out.file, sep = ' ')
  if(is.minMax) cmd <- paste(cmd, '-minMax', sep = ' ')

  # Run bigWigAverageOverBed
  print(cmd)
  system(cmd)

  # Read result file
  out.columns <- c('name', 'size', 'covered', 'sum',
                   'mean0', 'mean', 'min', 'max')
  bw.average <- read_tsv(out.file, col_names = out.columns) %>%
    arrange(name)

  return(bw.average)
}

#' Computes a metaprofile given GRange object and *.bw file
#'
#' @param region.gr
#' @param bw.file
#' @param region.names
#'
#' @return
#'
#' @import GenomicRanges
#' @importFrom rtracklayer import.bw
#' @importFrom dplyr data_frame %>% mutate group_by summarise
#' @importFrom tidyr unnest
#'
#' @export
#'
computeProfile <- function(region.gr, bw.file, region.names = NULL) {
  # Load ATAC-seq signal from bigwig
  atac.signal <- import.bw(bw.file, which = region.gr)
  # Compute overlap
  ov <- findOverlaps(atac.signal, region.gr)
  # Compute mean signal and std error
  atac.profiles <- data_frame(peak_id = subjectHits(ov),
                              peak_mid = mid(ranges(region.gr))[subjectHits(ov)],
                              bin_start = start(atac.signal)[queryHits(ov)],
                              bin_end = end(atac.signal)[queryHits(ov)],
                              bin_width = width(atac.signal)[queryHits(ov)],
                              score = atac.signal$score[queryHits(ov)],
                              name = region.names) %>%
    mutate(rel_start = bin_start - peak_mid,
           rel_end = bin_end - peak_mid) %>%
    dplyr::select(peak_id, name, rel_start, rel_end, score) %>%
    group_by(peak_id, name, rel_start, rel_end) %>%
    mutate(rel_pos = list(rel_start:rel_end)) %>%
    unnest() %>%
    group_by(name, rel_pos) %>%
    summarise(mean_score = mean(score),
              std_error = atacR:::stdError(score))
  return(atac.profiles)
}
