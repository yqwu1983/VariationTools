#!/usr/bin/env Rscript

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

# Given a GTF file of genes and a VCF file of variants, the script calculates
# the number of variants within the gene exons for each gene.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(VariantAnnotation))
source('~/Programs/VariationTools/track_routines.R')

parse.command.line.args <- function() {
  # Parses command-line arguments using routines from the argparse module.
  #
  # Parameters:
  #   None.
  #
  # Returns:
  #   A named list of script parameters.
  
  parser <- ArgumentParser()  
  parser$add_argument('gtf_file', type='character', 
                      help='a GTF file of annotated genes')
  parser$add_argument('vcf_file', type='character', 
                      help='a VCF file of genomic variants variants')
  parser$add_argument('genome_length', type='character',
                      help='a name of a file containing sequence lengths')
  parser$add_argument('--coding_transcripts', type='character',
                      help='a list of transcripts to be examined')
  parser$add_argument('--gene_names', type='character',
                      help='a file containing transcript and gene names')
  parser$add_argument('output', type='character', help='an output file')
  args <- parser$parse_args()
}

read.gene.transcripts <- function(file) {
  # Given a name of a file containing a table of gene and transcript IDs, read
  # the IDs from it and return a named vector which names are transcript IDs
  # and names are gene IDs.
  #
  # Args:
  #   file: a name of a file containing a table of gene and transcript IDs
  #
  # Return:
  #   A named vector of gene IDs which names are transcript IDs; it can be used
  #   to assign the corresponding gene IDs to the transcripts.
  
  gene.transcripts <- read.table(file, as.is=TRUE)
  names(gene.transcripts) <- c('TRANSCRIPT.ID', 'GENE.ID')
  
  result <- gene.transcripts$GENE.ID
  names(result) <- gene.transcripts$TRANSCRIPT.ID
  
  return(result)
}

# The main script routines ----

args <- parse.command.line.args()

message(paste('Reading GTF file', basename(args$gtf_file), '...'))
transcripts <- read.gtf(args$gtf_file, sequence.len.file=args$genome_length)
message('Reading and filtering protein-coding genes...')
if (!is.null(args$coding_transcripts)) {
  # remove transcripts of non-coding genes, if required
  coding.transcripts <- scan(args$coding_transcripts, what=character())
  transcripts <- transcripts[transcripts$transcript_id %in% coding.transcripts]
}
message('Reading and applying gene names...')
if (!is.null(args$gene_names)) {
  # assign gene names, if required
  transcript.genes <- read.gene.transcripts(args$gene_names)
  transcripts$gene_id <- transcript.genes[transcripts$gene_id]
}

# Now we split the transcripts by their gene names and reduce exon ranges for
# each gene.
gene.exons <- reduce(split(transcripts, transcripts$gene_id))

# Next, we read variants and for each gene, get the number of them within its
# exons.
message(paste('Reading variants from', basename(args$vcf_file), '...'))
seq.length <- read.seq.length(args$genome_length)
variants <- readVcfAsVRanges(args$vcf_file, 
                             Seqinfo(seq.length$NAME,
                                     seqlengths=seq.length$LENGTH))

message('Calculating exon SNV counts...')
variant.counts <- countOverlaps(gene.exons, variants)

message('Writing output...')
# Write the obtained counts to the specified output file.
write.table(variant.counts, args$output, row.names=TRUE, col.names=FALSE,
            quote=FALSE, sep='\t')
message('completed.')
