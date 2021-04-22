#!/usr/bin/env Rscript

## Script to read the Local PFAM Database output table into a CSV
## Note: Local PFAM Db output is different fromt he PFAM REST API output.

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_file"),
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
    default = NA,
    type = 'character',
    help = 'Output table containing json info'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## Functions ##
read.PFAM.domtblout <- function(file.path, transcript.id){
  ## Function to read the output of the Local installation of the PFAM database when --domtblout is specifief as output
  
  # 1. Names to add to the PFAM alignment table
  table.names <- c("transcript.id", "pfam.name", "pfam.accesion", "hmm.length", "transctipt.id", "accession.2", "query.seq.length", "E.value", "bit.score.overall", "bias.overall", "domain.number", "n.domain", "c.Evalue", "i.Evalue", "bit.score.domain", "bias.domain", "hmm.from", "hmm.to", "alignment.from", "alignment.to", "envelope.from", "envelope.to", "accuracy", "pfam.description", "pfam.match")
  # 2. Read 
  line.list = readLines(file.path)
  uncommented.lines <- grep("#", line.list)
  
  ## 2.1. PFAM Empty Alignment 
  if(length(uncommented.lines) == length(line.list)){
    empty.align <- rep(c(NA), times = length(table.names))
    names(empty.align) <- table.names
    ## Add transcript.id
    empty.align[["transcript.id"]] = transcript.id
    ## Add PFAM match logical value
    pfam.match = FALSE
    empty.align[["pfam.match"]] = pfam.match
    pfam.df = data.frame(t(empty.align))
    
  } else {
    
    ## 2.2. Alignment with content.
    all.lines <- 1:length(line.list)
    filled.lines <- all.lines[!all.lines %in% uncommented.lines]
    pfam.align.df = read.table(file.path, comment.char = "#", sep = "")
    ## Edits: last column ('pfam.description') is separated since it contains >1 words
    rm.desc.pfam.align.df <- pfam.align.df[, 1:22]
    pfam.desc.vec <- apply(pfam.align.df[, 23:ncol(pfam.align.df)], 1, paste , collapse = "_" )
    ### Concatenate back
    pfam.desc.df <- rm.desc.pfam.align.df
    pfam.desc.df[, 23] <- pfam.desc.vec
    rm(rm.desc.pfam.align.df)
    ## Add transcript.id
    transcript.id.df <- data.frame(Transcript.id = rep(transcript.id, times = nrow(pfam.align.df)))
    ## Add PFAM match logical value
    pfam.match = TRUE
    pfam.match.df <- data.frame(pfam.match = rep(pfam.match, times = nrow(pfam.align.df)))
    ## Merge dfs
    pfam.df <- cbind(transcript.id.df, pfam.desc.df, pfam.match.df)
    names(pfam.df) = table.names
  }
  # Close connection 
  return(pfam.df)
}


suppressPackageStartupMessages(require(rjson))

## 0. Params


# 1. Coerce dmtblout to table
pfam.table <- read.PFAM.domtblout(file.path = opt$input_file, transcript.id = opt$transcript_id)
# 2. Save table
write.table(pfam.table, file = opt$output_table,  quote = F, row.names = F, sep = "\t") # Important to save with a separating character that won't be present in the character strings output of the alignment, if not htis incurs problems in  the following Merging step
