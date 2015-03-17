#!/usr/bin/env Rscript

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

# The script processes output of compare_allele_counts.py script and adds the
# column of Fisher's exact test p-values to it.

suppressPackageStartupMessages(library(argparse))

parse.command.line.args <- function() {
  # Parses command-line arguments using routines from the argparse module.
  #
  # Parameters:
  #   None.
  #
  # Returns:
  #   A named list of script parameters.
  
  parser <- ArgumentParser()  
  parser$add_argument('input', type='character', help='an input file')
  parser$add_argument('output', type='character', help='an output file')
  args <- parser$parse_args()
}

add.fisher.pvalues <- function(df) {
  # Given a data frame from a file produced by compare_allele_counts.py script,
  # add a column of p-values for Fisher's exact test.
  #
  # Arguments:
  #   df: a data frame read from an output file of compare_allele_counts.py
  #
  # Returns:
  #   The column of Fisher's exact test p-values corresponding to each row of
  #   the specified data frame.
  
  names(df) <- c('CHROM', 'POS', 'REF', 'ALT', 'REF_A', 'ALT_A', 'REF_B',
                 'ALT_B')
  num.variants <- nrow(df)
  
  p.values <- numeric(num.variants)
  for (i in 1:num.variants) {
    contingency.table <- matrix(c(df$REF_A[i], df$ALT_A[i], 
                                  df$REF_B[i], df$ALT_B[i]),
                                byrow=TRUE, nrow=2, ncol=2)
    p.values[i] <- fisher.test(contingency.table)$p.value
  }
  
  return(p.values)
}

# The main script routines ----

args <- parse.command.line.args()
df <- read.table(args$input, as.is=TRUE)
df$P.VALUE <- add.fisher.pvalues(df)
write.table(df, args$out, col.names=FALSE, row.names=FALSE, quote=FALSE,
            sep='\t')
