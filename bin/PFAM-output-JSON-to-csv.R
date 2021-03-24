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
json.to.table <- function(json.file, transcript.id){
  ## Build df with extracted information from the PFAM alingment JSON file
  
  if(length(json.file) > 0){
    json.file = json.file[[1]]
    # edit some params
    if(json.file$sig == 1){ signif <- TRUE } else{  signif <- FALSE}
    if(is.null(json.file$act_site)){active_site <- NA}
    alignment.vec <- c(
      # Transcript ID
      transcript_id = transcript.id, 
      # Alignment
      pfam.name = json.file$name, pfam.accesion = json.file$acc,  pfam.description = json.file$desc, pfam.clan = json.file$clan,
      # Type
      pfam.type = json.file$type,
      # Signigicance
      significant = signif,
      # Envelope
      env.start = as.numeric(json.file$env$from), env.end = as.numeric(json.file$env$to), 
      # Alignment
      e.value = as.numeric(json.file$evalue), alignment.from =  as.numeric(json.file$seq$from), alignment.to =  as.numeric(json.file$seq$to),
      # HMM
      hmm.length = as.numeric(json.file$model_length), hmm.from = as.numeric(json.file$hmm$from),  hmm.to = as.numeric(json.file$hmm$to),
      # other data
      bit.score = as.numeric(json.file$bits), active.site = active_site, seq.name = json.file$seq$name, 
      pfam.match = TRUE
  )
  } else {
    alignment.vec <- c(
      transcript_id = transcript.id, 
      # alignment
      pfam.name = NA, pfam.accesion = NA,  pfam.description = NA, pfam.clan = NA, pfam.type = NA, significant = NA,
      env.start = NA, env.end = NA,  e.value = NA, alignment.from =  NA, alignment.to =  NA, hmm.length = NA, hmm.from = NA,
      hmm.to = NA, bit.score = NA, active.site = NA, seq.name = NA,  pfam.match = FALSE ) 
  }

  # build data.frame
  data.frame(as.list(alignment.vec)
  )
}


suppressPackageStartupMessages(require(rjson))

json.file = fromJSON(file = opt$input_json)
transcript.id = opt$transcript_id
# coerce json to table
pfam.table <- json.to.table(json.file = json.file, transcript.id = transcript.id)
# Save table
write.csv(pfam.table, file = opt$output_table,  quote = F, row.names = F)
