# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

# This script is a module containing routines to perform SNV density analysis.

suppressPackageStartupMessages(library(GenomicRanges))

get.windows <- function(regions, window.len,
                        remainder.threshold=0.5*window.len, verbose=FALSE) {
  # Given a GRanges object representing regions in a genome, return the GRanges
  # object of genome windows of the specified length distributed in the
  # regions.
  #
  # Arguments:
  #   regions: a GRanges object of regions in a genome
  #   window.len: a length of windows to be contructed for the specified genome
  #     regions
  #   remainder.threshold: the least length of a window situated in the end of
  #     a region to be included in the result
  #   verbose: show a progress bar indicating the window computation progress
  #
  # Returns:
  #   The GRanges object of the constructed windows within the regions.
  
  result <- GRanges(seqinfo=seqinfo(regions))
  
  if (verbose) {
    message(paste('Calculating windows of', window.len, 'bp in', 
                  length(regions), 'regions...'))
    pb <- txtProgressBar(min=0, max=length(regions), style=3)
  }
  
  # remove the regions which length is less than the window size
  regions <- regions[width(regions) >= window.len]
  
  # We create vectors of sequence names, start and end positions for the given
  # ranges to accelerate the loop below.
  region.seq.names <- seqnames(regions)
  region.start <- start(regions)
  region.end <- end(regions)
  for (i in 1:length(regions)) {
    temp.start <- seq(region.start[i], region.end[i], by=window.len)
    temp <- GRanges(seqnames=region.seq.names[i],
                    ranges=IRanges(start=temp.start, width=window.len),
                    seqinfo=seqinfo(regions))
    # check the last window of the region
    if (width(temp[length(temp)]) < remainder.threshold) {
      temp <- temp[-length(temp)]
    }
    result <- c(result, temp)
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) {
    close(pb)
  }
  
  return(result)
}

get.feature.counts <- function(regions, features) {
  # Given two GRanges objects of genomic regions and features, return the 
  # GRanges object containing the additional column of feature counts within
  # the regions.
  #
  # Arguments:
  #   regions: a GRanges object of regions in a genome
  #   features: a GRanges object of genome features
  #
  # Returns:
  #   The GRanges object of genome regions with the column of feature counts
  #   within them.
  
  counts <- countOverlaps(regions, features, type="within")
  
  mcols(regions) <- data.frame(feature.counts=counts)
  
  return(regions)
}
