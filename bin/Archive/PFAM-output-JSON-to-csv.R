#!/usr/bin/env Rscript

## Script to read the PFAM output json fail into an R table

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_json"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Json object output by PFAM'
  ),
  make_option(
    c("-t", "--transcript_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Transcript ID of the sequence that was aligned to PFAM.'
  ),
  make_option(
    c("-o", "--output_table"),
    action = "store",
    default = "",
    type = 'character',
    help = 'Output table containing json info'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## Functions ##
build.pfam.alignment.df <- function(pfam.hit, transcript.id){
  # Build df from a single PFAM alignment 
  
  # edit some params
  if(pfam.hit$sig == 1){ signif <- TRUE } else{  signif <- FALSE}
  if(is.null(pfam.hit$act_site)){active_site <- NA}
  
  alignment.vec <- c(
    # Transcript ID
    transcript_id = transcript.id,
    # Alignment
    pfam.name = pfam.hit$name, pfam.accesion = pfam.hit$acc,  pfam.description = I(pfam.hit$desc), pfam.clan = pfam.hit$clan,
    # Type
    pfam.type = pfam.hit$type,
    # Significance
    significant = signif,
    # Envelope
    env.start = as.numeric(pfam.hit$env$from), env.end = as.numeric(pfam.hit$env$to),
    # Alignment
    e.value = as.numeric(pfam.hit$evalue), alignment.from =  as.numeric(pfam.hit$seq$from), alignment.to =  as.numeric(pfam.hit$seq$to),
    # HMM
    hmm.length = as.numeric(pfam.hit$model_length), hmm.from = as.numeric(pfam.hit$hmm$from),  hmm.to = as.numeric(pfam.hit$hmm$to),
    bit.score = as.numeric(pfam.hit$bits), active.site = active_site, seq.name = pfam.hit$seq$name,
    pfam.match = TRUE)
    # Return horizontal df
    return(data.frame(t(alignment.vec)))
}

json.to.table <- function(json.file, transcript.id){
  ## Build a df summarizing the output of a PFAM alignment json file
  ## PFAM alignment JSON file can contain: one, multiple or none hits.
  
  # PFAM Hit or Multi-Hit
  if(length(json.file) > 0){
    hit.df.list = lapply(json.file, function(pfam.hit) build.pfam.alignment.df(pfam.hit, transcript.id=transcript.id))
    # Merge df list to DF
    pfam.df = Reduce(rbind, hit.df.list)
    
    # PFAM Empty hit
  } else {
    empty.vec <- c(
      transcript_id = transcript.id,
      # alignment
      pfam.name = NA, pfam.accesion = NA, pfam.description = NA, pfam.clan = NA, pfam.type = NA, significant = NA,
      env.start = NA, env.end = NA, e.value = NA, alignment.from = NA, alignment.to = NA, hmm.length = NA, hmm.from = NA,
      hmm.to = NA, bit.score = NA, active.site = NA, seq.name = NA, pfam.match = FALSE )
    # build data.frame
    pfam.df = data.frame(as.list(empty.vec))
  }
  return(noquote(pfam.df))
}


suppressPackageStartupMessages(require(rjson))

json.file = fromJSON(file = opt$input_json)
transcript.id = opt$transcript_id
# coerce json to table
pfam.table <- json.to.table(json.file = json.file, transcript.id = transcript.id)
# Save table
write.table(pfam.table, file = opt$output_table,  quote = F, row.names = F, sep = "\t") # Important to save with a separating character that won't be present in the character strings output of the alignment, if not htis incurs problems in  the following Merging step
