#!/usr/bin/env Rscript

## Hmmscan Domtblout Parser
## This script utilises functions from the rhmmer CRAN package [https://github.com/arendsee/rhmmer/issues]

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--hmmscan_tbl"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Hmmscan output in --domtblout format'
  ),
  make_option(
    c("-t", "--transcript_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Transcript ID of the sequence that was aligned with HMMSCAN to PFAM database'
  ),
  make_option(
    c("-r", "--rhmmer_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the "rhmmer.R" script (rhmmer CRAN package), from which we are using hmmscan output parsing functions.'
  ),
  make_option(
    c("-o", "--output_table"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output table containing json info'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))


suppressPackageStartupMessages(require(magrittr))

# 2) Coerce tbl to df
tbl.to.df <- function(tb, transcript_id){
  ## Coerce hmmscan parsed tbl (from 'read_domtblout()') to df. Replenish with NAs if tb is empty
  # If empty tb
  if(nrow(tb) == 0){
    # create empty df
    df <- as.data.frame(setNames(rep(list(NA), length(names(tb))), names(tb)))
    # Add transcript name
    df[["query_name"]] <- transcript_id
    df[["pfam_match"]] <- FALSE
  } else {
    df <- data.frame(tb)
    df[["pfam_match"]] <- TRUE
  }
  return(df)
}

# 3)  [NOT RUNNING FOR NOW - NOT CLEAR, Some fractions are above 1, Potential wrong formula] Calculate domain coverage 
## In addition, this can be included in a data analysis step, it is out of the scope of a simple 'parser'
frac.alig.domain <- function(df){
  ## Calculate the fraction of AA of the HMM model spanned by the PFAM alignment
  ## This is to address Partiality in PFAM alignments. Based on the idea that an alignment of the query sequence (or a part of the query sequence) against 100% of the PFAM domain can be considered a functional domain. Whereas an alignment of 30% of the domain, can't.
  df[["partiality"]] = (df[["hmm_to"]] + 1 - df[["hmm_from"]])/df[["domain_len"]]
  return(df)
}

## RUN ## 
# 0. rhmmer (CRAN) functions
source(opt$rhmmer_path)
# 1. Parser HMMScan Output Table
hmmscan.tb <- read_domtblout(opt$hmmscan_tbl)
# 2. Coerce to df
hmmscan.df <- tbl.to.df(hmmscan.tb, transcript_id = opt$transcript_id)
# 3. Add fraction of alignment
hmmscan.df <- frac.alig.domain(hmmscan.df)
# 4. Save TSV
write.table(hmmscan.df, file = opt$output_table, quote = F, sep = "\t")
