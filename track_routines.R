# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

# This script is a module containing routines to handle genome tracks.

suppressPackageStartupMessages(library(GenomicRanges))

read.sequence.len <- function(file) {
  # Read lengths of sequences from the specified file and return them as a data
  # frame of two columns: sequence names and their lengths in base pairs.
  #
  # Args:
  #   file: a name of a file containing sequence lengths
  #
  # Returns:
  #   The data frame of two columns (NAME and LENGTH) that contains names of
  #   genomic sequences and their lengths in base pairs.
  #
  # Note:
  #   The specified file must contain two columns: sequence names and lengths. 
  #   One may obtain it for sequences in a FASTA file using faSize from Jim
  #   Kent's utilites with the -detailed option.
  
  result <- read.table(file, as.is=TRUE)
  if (length(result) != 2) {
    stop('the incorrect number of columns in a sequence length file')
  }
  names(result) <- c('NAME', 'LENGTH')
  if (!is.numeric(result$LENGTH)) {
    stop('incorrect sequence length values')
  }
  
  return(result)
}

read.bed <- function(file, only.loci=FALSE, sequence.len.file=NULL) {
  # Read contents of a BED track file and return the GRanges object
  # representing them.
  #
  # Args:
  #   file: a name of a track file in the BED format
  #   only.loci: get only positions of features in a genome without any
  #     additional information
  #   sequence.len.file: a name of a file containing lengths of genome
  #     sequences of the track
  #
  # Returns:
  #   The GRanges object containing genomic positions from the specified file.
  #
  # Note:
  #   The name of a sequence length file is an optinal argument. However, it is
  #   recommended to use it if possible to ensure that regions in a track file
  #   are correct.
  
  bed.names <- c('CHROM',
                 'CHROM.START',
                 'CHROM.END',
                 'NAME',
                 'SCORE',
                 'STRAND',
                 'THICK.START',
                 'THICK.END',
                 'ITEM.RGB',
                 'BLOCK.COUNT',
                 'BLOCK.SIZES',
                 'BLOCK.STARTS')
  
  track <- read.table(file, as.is=TRUE)
  names(track)[1:min(length(track), length(bed.names))] <- 
    bed.names[1:min(length(track), length(bed.names))]
  
  track.seqinfo <- NULL
  if (!is.null(sequence.len.file)) {
    seq.len <- read.sequence.len(sequence.len.file)
    track.seqinfo <- Seqinfo(seq.len$NAME, seqlengths=seq.len$LENGTH)
  } else {
    warning('a sequence length file is not specified')
  }
  
  # we convert positions to be 1-based and inclusive
  track$CHROM.START <- track$CHROM.START + 1
  
  result <- GRanges(seqnames=track$CHROM,
                    ranges=IRanges(start=track$CHROM.START,
                                   end=track$CHROM.END),
                    seqinfo=track.seqinfo)
  
  # if the strand column is specified, use it and remove from the data frame
  if ('STRAND' %in% names(track)) {
    strand(result) <- track$STRAND
    track$STRAND <- NULL
  }
  
  # process optional BED fields
  if (!only.loci && length(track) > 3) {
    mcols(result) <- track[, 4:length(track)]
  }
  
  return(result)
}
