#!/usr/bin/env Rscript

## This script filters the PFAM output as defined by our filters

# Vulnerabilities: We are generating 2 objects (when we could be using one), this takes more memory.

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--pfam_alignment"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'PFAM alignment file'
  ), 
  # Thresholds
  make_option(
    c("-s", "--sequence_evalue"),
    action = "store",
    default = 10,
    type = 'numeric',
    help = 'Sequence E-value threshold. The E-value is the number of hits that would be expected to have a score equal to or better than this value by chance alone. A good E-value is much less than 1.'
  ), 
  make_option(
    c("-c", "--sequence_score"),
    action = "store",
    default = 1e-03,
    type = 'numeric',
    help = 'Sequence Bit Score. Represents the likelihood of alignment being emitted by the model.'
  ), 
  make_option(
    c("-d", "--domain_evalue"),
    action = "store",
    default = 1e-03,
    type = 'numeric',
    help = 'Domain E-value threshold. This value represents both the conditional-E-Value (smaller seach space, more permissive significance measure) and the independent-E-value (on entire database, more stringent significant measure)'
  ), 
  make_option(
    c("-e", "--domain_score"),
    action = "store",
    default = 10,
    type = 'numeric',
    help = 'Domain Bit Score. Represents the likelihood of alignment being emitted by the model.'
  ), 
  make_option(
    c("-a", "--accuracy"),
    action = "store",
    default = 0.8,
    type = 'numeric',
    help = 'The mean posterior probability of aligned residues in the MEA (maximum expected accuracy) alignment; a measure of how reliable the overall alignment is (from 0 to 1, with 1.00 indicating a completely reliable alignment
according to the model).'
  ), 
  make_option(
    c("-p", "--partiality"),
    action = "store",
    default = 0.9,
    type = 'numeric',
    help = 'Segment of the domain covered by the alignment. We consider an a alignment a "high-confidence domain" if it spans 90% of the domain length.'
  ), 
  #
  make_option(
    c("-o", "--pfam_filtered"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Filtered PFAM Alignment according to our defined filters.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# 0. Params 
sequence_evalue <- opt$sequence_evalue
sequence_score <- opt$sequence_score
domain_evalue <- opt$domain_evalue
domain_score <- opt$domain_score
accuracy <- opt$accuracy
partiality <- opt$partiality

# 0. Checks: Are all values numeric ? 
threshold.vec <- setNames(c(sequence_evalue, sequence_score, domain_evalue, domain_score, accuracy, partiality),
                          nm = c("sequence_evalue", "sequence_score", "domain_evalue", "domain_score", "accuracy", "partiality"))
for(i in 1:length(threshold.vec)){
  if(!is.numeric(threshold.vec[i])){
    print(paste0("Threshold:\n", threshold.vec[i], "\n is not numeric. Revise!"))
  }
}

## 1. Read file
pfam.df <- read.table(opt$pfam_alignment, header = T, sep = "\t", quote = "\"")
# Get transcript
transcript.id <- unique(pfam.df[["query_name"]])

# 2. Filter
if(all(pfam.df[["pfam_match"]])){
  
  pfam.df.filt <- pfam.df[pfam.df[["sequence_evalue"]] <= threshold.vec[["sequence_evalue"]] & # Sequence E-Value
            pfam.df[["sequence_score"]] >= threshold.vec[["sequence_score"]] & # Sequence Bit-score
            pfam.df[["domain_cevalue"]] <= threshold.vec[["domain_evalue"]] & # Domain conditional E-value
            pfam.df[["domain_ievalue"]] <= threshold.vec[["domain_evalue"]] & # Domain independent E-value
            pfam.df[["domain_score"]] >= threshold.vec[["domain_score"]] & # Domain Bit-score
            pfam.df[["acc"]] >= threshold.vec[["accuracy"]] & # Alignment accuracy
            pfam.df[["partiality"]] >= threshold.vec[["partiality"]]  # Partiality
            , ]
  ## If zero-row df, fill with NA
  if(nrow(pfam.df.filt) == 0){
    pfam.df.filt[1, ] <- NA
    pfam.df.filt[["query_name"]] <- transcript.id
    pfam.df.filt[["pfam_match"]] <- FALSE
    }

} else if (all(pfam.df[["pfam_match"]]) == F) { 
  pfam.df.filt <- pfam.df
}

# 3. Save
write.table(pfam.df.filt, file = opt$pfam_filtered, sep = "\t", quote = F, col.names = T, row.names = F)