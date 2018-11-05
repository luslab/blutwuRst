#' Load interactions from *.bed files and save as RData file
#'
#' @param interA.file
#' @param interB.file
#' @param inter.meta
#'
#' @import GenomicInteractions
#'
#' @return
#'
prepareInteractionDb <- function(interA.file, interB.file, inter.meta, out.rdata) {
  # Load interaction database metadata
  interDb.meta <- read_tsv(inter.meta)

  # Load interactions from *.bed file
  interA.gr <- import.bed(interA.file)
  interB.gr <- import.bed(interB.file)

  # Combine interactions and metadata in GenomicInteractions object
  interactions <- GenomicInteractions(interA.gr, interB.gr)
  mcols(interactions)  <- interDb.meta

  # Save interactions object to RData file
  save(interactions, file = out.rdata)
}

#' Load interaction database (4DGenome)
#'
#' @param genome
#'
#' @return
#'
#' @import GenomicInteractions
#'
#' @export
#'
loadInteractionDb <- function(genome = 'mm10') {
  # Find
  if (genome == 'mm10') {
    rdata.file <- system.file('extdata', '4dgenome_interactions_mm10.RData',
                             package = 'blutwuRst')
  } else{
    print(sprintf('No data available for %s', genome))
  }
  # Load interaction data from RData file
  load(rdata.file)
  return(interactions)
}

#' Subsets interactions by a given regions
#' (--> get interactions within region)
#'
#' @param interactions
#' @param region
#'
#' @return
#'
#' @import GenomicInteractions
#'
#' @export
#'
subsetInteractionsByRegion <- function(interactions, region) {
  ov1 <- findOverlaps(anchorOne(interactions), region) %>% queryHits(.)
  ov2 <- findOverlaps(anchorTwo(interactions), region) %>% queryHits(.)
  i.ov <- c(ov1, ov2) %>% table(.)
  i.ov <- names(i.ov)[i.ov > 1] %>% as.numeric(.)
  return(interactions[i.ov,])
}

#' Subset interactions by a given gene target
#'
#' @param interactions
#' @param gene
#'
#' @return
#'
#' @import GenomicInteractions
#'
#' @export
#'
subsetInteractionsByGene <- function(interactions, gene) {
  i.gene <- c(grep(interactions$Agene, pattern = gene),
              grep(interactions$Bgene, pattern = gene))
  return(interactions[i.gene,])
}


