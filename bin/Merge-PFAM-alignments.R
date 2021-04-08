#!/usr/bin/env Rscript

# Pre process Alternative Splicing (AS) input files
## Output various things: i) df with functional concordance, ii) df with no-functional concordance, iii) list of unique transcript ids for PFAM query

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_pfam"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path directory containing all the PFAM alignment csv files.'
  ),
  make_option(
    c("-o", "--output_pfam_df"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Merged dataframe with all alignments from PFAM.'
  )
)


opt <- parse_args(OptionParser(option_list=option_list))

# 1. List files in dir
files = list.files(opt$input_pfam, pattern = ".txt", full.names = T)
file.list = lapply(files, read.table, sep = "\t", header = T)
# 2. Merge files to single df
pfam.df = Reduce(rbind, file.list)
# 3. Save file
write.table(pfam.df, file = opt$output_pfam_df, sep = "\t", quote = F, row.names = T)
