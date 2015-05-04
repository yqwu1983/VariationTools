#!/bin/env Rscript

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

# This script draws ideogram representing homozygosity regions for the
# specified genome from the BED file of SNV frequencies in genome windows. It
# implents the approach described in
#
# Pontius, Joan U., James C. Mullikin, Douglas R. Smith, Kerstin Lindblad-Toh, 
# Sante Gnerre, Michele Clamp, Jean Chang et al. "Initial sequence and
# comparative analysis of the cat genome." Genome research 17, no. 11 (2007):
# 1675-1689.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Gviz))

parse.command.line.args <- function() {
  # Parses command-line arguments using routines from the argparse module.
  #
  # Parameters:
  #   None.
  #
  # Returns:
  #   A named list of script parameters.
  
  parser <- ArgumentParser()
  parser$add_argument('snv_freq_bed',
                      help='a BED file of SNV frequencies in genome windows')
  parser$add_argument('genome', help='a genome which variants are analyzed')
  parser$add_argument('gaps', help='a BED file of gaps in the genome')
  parser$add_argument('output_dir', 
                      help='an output directory for ideogram files')
  parser$add_argument('-c', '--count', type='integer', default=2,
                      help='the maximum number of SNVs in a window to consider it homozygous')
  args <- parser$parse_args()
}

read.snv.counts <- function(filename) {
  # Read SNV counts in genome regions from the specified file.
  #
  # Args:
  #   filename: a name of a BED3 file with the fourth column of SNV counts
  #
  # Return:
  #   A data frame of four columns: a sequence, start and end positions and the
  #   number of SNVs in the sequence region.
  
  temp <- read.table(filename, as.is=TRUE)
  names(temp) <- c('SEQ', 'START', 'END', 'SNV.COUNT')
  temp$START <- temp$START + 1
  return(temp)
}

get.gaps <- function(filename) {
  # Given a name of a BED file with gap regions, read them and return as a
  # GRanges object.
  #
  # Args:
  #   filename: a name of a BED file with gap regions
  #
  # Return:
  #   A GRanges object representing gap regions from the specified BED file.
  
  temp <- read.table(filename, as.is=TRUE)
  names(temp) <- c('SEQ', 'START', 'END')
  temp$START <- temp$START + 1
  result <- GRanges(seqnames=temp$SEQ, 
                    ranges=IRanges(start=temp$START, end=temp$END))
  return(result)
}

identify.homozygosity.regions <- function(df, gap.regions, max.snvs) {
  # Given a data frame of SNV counts within genome regions and the maximum
  # number of SNVs for a region to be considered homozygous, return a GRanges
  # object with continious ranges of homozygosity or non-homozygosity regions
  # with additional column indicating the region type.
  #
  # Arguments:
  #   df: a data frame with four columns representing SNV counts in genome
  #     regions, like the one produced by the read.snv.counts function
  #   gap.regions: a GRanges object of gap regions
  #   max.snvs: the maximum number of SNVs in a region to be considered
  #     homozygous
  #
  # Return:
  #   A GRanges object of reduced ranges of the genome regions with additional
  #   column indicating whether the region is homozygous or not.
  
  df$IS.HMZ <- (df$SNV.COUNT <= max.snvs)
  df <- GRanges(seqnames=df$SEQ, ranges=IRanges(start=df$START, end=df$END),
                is.homozygous=df$IS.HMZ)
  if (!is.null(gap.regions)) {
    df <- df[!as.logical(countOverlaps(df, gap.regions, type="any"))]    
  }
  df <- split(df, df$is.homozygous)
  df <- reduce(df)
  result <- list()
  for (i in c("TRUE", "FALSE")) {
    result[[i]] <- GRanges(seqnames=seqnames(df[[i]]), 
                           ranges=IRanges(start=start(df[[i]]), 
                                          end=end(df[[i]])),
                           is.homozygous=rep(as.logical(i), length(df[[i]])))
  }
  # concatenate the regions
  result <- do.call(c, unname(as.list(result)))
  # sort the regions
  result <- sort(result)
  return(result)
}

plot.hmz.regions <- function(genome.regions, genome, output.dir, 
                             width=200, height=100) {
  # Given a GRanges object with the additional field indicating if the regions
  # are homozygous or not, produce the ideogram plot.
  #
  # Args:
  #   genome.regions: a GRanges object indicating homozygosity regions
  #   genome: a name of a genome which chromosomes are considred
  #   output.dir: an output directory for the figures
  #   width: width of the figure in pixels
  #   height: height of the figure in pixels
  #
  # Returns:
  #   None.
  
  genome.regions <- split(genome.regions, seqnames(genome.regions))
  chromosome.names <- grep('_', names(genome.regions), invert=TRUE, value=TRUE)
  max.chr.size <- max(unlist(lapply(genome.regions, 
                                    function(x) end(x)[length(x)])))
  for (i in chromosome.names) {
    chr.size <- end(genome.regions[[i]])[length(genome.regions[[i]])]
    png(file.path(output.dir, paste(i, '.png', sep='')), 
        round(width*chr.size/max.chr.size), 10, units="mm", res=300)
    track <- IdeogramTrack(genome=genome, chromosome=i)
    pb <- txtProgressBar(min=0, max=length(genome.regions[[i]]), style=3)
    for (j in 1:length(genome.regions[[i]])) {
      cur.region <- genome.regions[[i]][j]
      plotTracks(track, from=start(cur.region),
                 to=end(cur.region),
                 fill=ifelse(cur.region$is.homozygous, 'green', 'red'),
                 col=ifelse(cur.region$is.homozygous, 'green', 'red'), 
                 add=(i > 1))
      setTxtProgressBar(pb, j)
    }
    close(pb)
    dev.off()
  }
}

args <- parse.command.line.args()

# args <- list(snv_freq_bed='~/Dropbox/Genome studies/felidae_genomes/data/2015-04-16/cheetah_hmz/L2_B1_chr_hmz.txt',
#              genome='felCat5',
#              gaps='~/Dropbox/Genome studies/felidae_genomes/data/2015-01-16/cat_gaps/cat_gaps.bed',
#              output_dir='~/Desktop/Cheetah_ideograms/'
#              count=40)

snv.counts <- read.snv.counts(args$snv_freq_bed)
if (!is.null(args$gaps)) {
  gap.regions <- get.gaps(args$gaps)
  gap.regions <- gap.regions[width(gap.regions) > 10^5]  
}
genome.regions <- identify.homozygosity.regions(snv.counts, gap.regions,
                                                args$count)
genome.regions <- genome.regions[grep('_', seqnames(genome.regions), 
                                      invert=TRUE)]
genome.regions <- genome.regions[seqnames(genome.regions) != 'chrM']
plot.hmz.regions(genome.regions, args$genome, args$output_dir)
